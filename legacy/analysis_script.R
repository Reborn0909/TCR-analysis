










setwd("DATA_ROOT")


dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)


packages_needed <- c("Seurat", "Signac", "tidyverse", "patchwork",
                     "ggplot2", "ggpubr", "future", "data.table")

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 8000 * 1024^2)


set.seed(42)

cat("，：", getwd(), "\n")






data_dir <- file.path(getwd(), "data")
setwd(data_dir)


zenodo_base <- "https://zenodo.org/records/7555405/files"


files_to_download <- c(
  "PBMC_vaccine_CITE.rds",
  "PBMC_vaccine_ASAP.rds",
  "PBMC_vaccine_ECCITE_TCR.rds",
  "antigen_module_genes.rds",
  "antigen_module_peaks.rds",
  "PBMC_vaccine_ASAP_fragments.tsv.gz",
  "PBMC_vaccine_ASAP_fragments.tsv.gz.tbi"
)


setwd("DATA_ROOT")





library(Seurat)
library(tidyverse)
library(patchwork)


cat("\n CITE-seq ...\n")
PBMC_CITE <- readRDS("data/PBMC_vaccine_CITE.rds")


cat("\n===== CITE-seq  =====\n")
cat("：", ncol(PBMC_CITE), "\n")
cat("：", nrow(PBMC_CITE[["RNA"]]), "\n")
cat("：", nrow(PBMC_CITE[["ADT"]]), "\n")


cat("\n===== celltypel2  =====\n")
print(table(PBMC_CITE$celltypel2))


cd8_types <- unique(PBMC_CITE$celltypel2)
cd8_types <- cd8_types[grepl("CD8", cd8_types, ignore.case = TRUE)]
cat("\n CD8⁺ T ：\n")
print(cd8_types)


cd8_cite <- subset(PBMC_CITE, subset = celltypel2 %in% cd8_types)

cat("\n CD8⁺ T ：", ncol(cd8_cite), "\n")


cat("\n===== CD8⁺ T （celltypel3） =====\n")
print(table(cd8_cite$celltypel3))


saveRDS(cd8_cite, "data/cd8_cite_subset.rds")


rm(PBMC_CITE)
gc()

cat("\n✅ CITE-seq CD8⁺ T \n")






cd8_cite <- readRDS("data/cd8_cite_subset.rds")

cat("CD8⁺ T ：", ncol(cd8_cite), "\n")





cat("\n...\n")

DefaultAssay(cd8_cite) <- "RNA"
cd8_cite <- NormalizeData(cd8_cite)
cd8_cite <- FindVariableFeatures(cd8_cite, nfeatures = 2000)
cd8_cite <- ScaleData(cd8_cite)
cd8_cite <- RunPCA(cd8_cite, npcs = 50)


ElbowPlot(cd8_cite, ndims = 50)
ggsave("results/figures/elbow_plot_CD8.pdf", width = 8, height = 6)


cd8_cite <- RunUMAP(cd8_cite, dims = 1:30)
cd8_cite <- FindNeighbors(cd8_cite, dims = 1:30)
cd8_cite <- FindClusters(cd8_cite, resolution = 0.6)

cat("：\n")
print(table(cd8_cite$seurat_clusters))


p1 <- DimPlot(cd8_cite, group.by = "seurat_clusters", label = TRUE, label.size = 5) +
  ggtitle("CD8⁺ T ")
ggsave("results/figures/umap_clusters_CD8.pdf", p1, width = 10, height = 8)

p2 <- DimPlot(cd8_cite, group.by = "celltypel3", label = TRUE, label.size = 4) +
  ggtitle("CD8⁺ T （）")
ggsave("results/figures/umap_celltypel3_CD8.pdf", p2, width = 12, height = 8)

p3 <- DimPlot(cd8_cite, group.by = "timepoint", label = TRUE) +
  ggtitle("")
ggsave("results/figures/umap_timepoint_CD8.pdf", p3, width = 10, height = 8)





DefaultAssay(cd8_cite) <- "ADT"
adt_features <- rownames(cd8_cite[["ADT"]])


cd27_name <- grep("^CD27$|CD27[^0-9]", adt_features, value = TRUE, ignore.case = TRUE)
cd27_clean <- cd27_name[!grepl("CD27[0-9]", cd27_name, ignore.case = TRUE)]
cd27_feature <- cd27_clean[1]
hla_dr_feature <- "HLA-DR"

cat(" HLA-DR ：", hla_dr_feature, "\n")
cat(" CD27 ：", cd27_feature, "\n")


hla_dr_expr <- FetchData(cd8_cite, vars = hla_dr_feature)[, 1]
cd27_expr <- FetchData(cd8_cite, vars = cd27_feature)[, 1]


hla_dr_threshold <- quantile(hla_dr_expr[hla_dr_expr > 0], 0.66, na.rm = TRUE)
cd27_threshold <- quantile(cd27_expr[cd27_expr > 0], 0.66, na.rm = TRUE)


cd8_cite$HLA_DR_status <- ifelse(hla_dr_expr > hla_dr_threshold, "HLA-DR+", "HLA-DR-")
cd8_cite$CD27_status <- ifelse(cd27_expr > cd27_threshold, "CD27+", "CD27-")
cd8_cite$HLA_DR_CD27_group <- paste(cd8_cite$HLA_DR_status, cd8_cite$CD27_status, sep="_")

cat("\n===== HLA-DR/CD27  =====\n")
print(table(cd8_cite$HLA_DR_CD27_group))

group_pct <- prop.table(table(cd8_cite$HLA_DR_CD27_group)) * 100
cat("\n：\n")
print(round(group_pct, 2))


p5 <- DimPlot(cd8_cite, group.by = "HLA_DR_CD27_group",
              cols = c("HLA-DR+_CD27+" = "#E41A1C",
                       "HLA-DR+_CD27-" = "#FF7F00",
                       "HLA-DR-_CD27+" = "#4DAF4A",
                       "HLA-DR-_CD27-" = "#999999")) +
  ggtitle("HLA-DR and CD27 expression groups")
ggsave("results/figures/umap_HLA_DR_CD27_groups.pdf", p5, width = 10, height = 8)

p6 <- FeaturePlot(cd8_cite, features = c(hla_dr_feature, cd27_feature),
                  blend = TRUE, blend.threshold = 0.5, ncol = 1)
ggsave("results/figures/featureplot_HLA_DR_CD27.pdf", p6, width = 8, height = 12)


saveRDS(cd8_cite, "data/cd8_cite_analyzed.rds")

cat("\n✅ \n")





library(ggpubr)


cd8_cite <- readRDS("data/cd8_cite_analyzed.rds")





antigen_genes <- readRDS("data/antigen_module_genes.rds")
cat("VI-GEM ", length(antigen_genes), "\n")

DefaultAssay(cd8_cite) <- "RNA"
cd8_cite <- AddModuleScore(cd8_cite,
                           features = list(VI_GEM = antigen_genes),
                           name = "VI_GEM")
colnames(cd8_cite@meta.data)[grepl("VI_GEM", colnames(cd8_cite@meta.data))] <- "VI_GEM_score"





DefaultAssay(cd8_cite) <- "RNA"

tf_genes <- c("FOXP3", "IKZF2", "TCF7", "BATF3", "EOMES", "TBX21")
co_stim_genes <- c("CD28", "ICOS", "TNFRSF9", "TNFRSF4")
co_inhib_genes <- c("PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2")

all_genes <- c(tf_genes, co_stim_genes, co_inhib_genes)
available_genes <- intersect(all_genes, rownames(cd8_cite[["RNA"]]))

cat("\n：", paste(available_genes, collapse = ", "), "\n")

for (gene in available_genes) {
  cd8_cite@meta.data[[gene]] <- FetchData(cd8_cite, vars = gene)[, 1]
}





comparison_data <- cd8_cite@meta.data %>%
  filter(HLA_DR_CD27_group %in% c("HLA-DR+_CD27+", "HLA-DR-_CD27-")) %>%
  select(HLA_DR_CD27_group, VI_GEM_score, all_of(available_genes))

cat("\n：\n")
print(table(comparison_data$HLA_DR_CD27_group))





plot_comparison <- function(data, feature, title = NULL) {
  data_sub <- data[!is.na(data[[feature]]), ]
  if (nrow(data_sub) == 0) return(NULL)

  ggplot(data_sub, aes(x = HLA_DR_CD27_group, y = .data[[feature]], fill = HLA_DR_CD27_group)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("HLA-DR+_CD27+" = "#E41A1C", "HLA-DR-_CD27-" = "#999999")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = feature, title = title) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5)
}





p_vigem <- plot_comparison(comparison_data, "VI_GEM_score", "VI-GEM module score")
ggsave("results/figures/comparison_VI_GEM.pdf", p_vigem, width = 6, height = 6)





tf_available <- intersect(tf_genes, available_genes)
p_tf_list <- list()
for (gene in tf_available) {
  p_tf_list[[gene]] <- plot_comparison(comparison_data, gene, gene)
}
p_tf_list <- p_tf_list[!sapply(p_tf_list, is.null)]
if (length(p_tf_list) > 0) {
  p_tf <- wrap_plots(p_tf_list, ncol = min(3, length(p_tf_list)))
  ggsave("results/figures/comparison_TF.pdf", p_tf, width = 15, height = 10)
}





co_inhib_available <- intersect(co_inhib_genes, available_genes)
p_coinhib_list <- list()
for (gene in co_inhib_available) {
  p_coinhib_list[[gene]] <- plot_comparison(comparison_data, gene, gene)
}
p_coinhib_list <- p_coinhib_list[!sapply(p_coinhib_list, is.null)]
if (length(p_coinhib_list) > 0) {
  p_coinhib <- wrap_plots(p_coinhib_list, ncol = min(3, length(p_coinhib_list)))
  ggsave("results/figures/comparison_co_inhib.pdf", p_coinhib, width = 15, height = 10)
}





stat_results <- data.frame(
  Feature = character(),
  Mean_DRpos_CD27pos = numeric(),
  Mean_DRneg_CD27neg = numeric(),
  Log2FC = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

features_to_test <- c("VI_GEM_score", tf_available, co_inhib_available)

for (feature in features_to_test) {
  if (feature %in% colnames(comparison_data)) {
    x <- comparison_data[[feature]][comparison_data$HLA_DR_CD27_group == "HLA-DR+_CD27+"]
    y <- comparison_data[[feature]][comparison_data$HLA_DR_CD27_group == "HLA-DR-_CD27-"]

    if (length(x) > 1 && length(y) > 1) {
      test <- wilcox.test(x, y)
      mean_x <- mean(x, na.rm = TRUE)
      mean_y <- mean(y, na.rm = TRUE)

      stat_results <- rbind(stat_results, data.frame(
        Feature = feature,
        Mean_DRpos_CD27pos = mean_x,
        Mean_DRneg_CD27neg = mean_y,
        Log2FC = log2((mean_x + 0.001) / (mean_y + 0.001)),
        P_value = test$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(stat_results) > 0) {
  stat_results$FDR <- p.adjust(stat_results$P_value, method = "BH")
  stat_results <- stat_results[order(stat_results$Log2FC, decreasing = TRUE), ]
}

cat("\n=====  =====\n")
print(stat_results)

write.csv(stat_results, "results/tables/stat_comparison_HLA_DR_CD27.csv", row.names = FALSE)





saveRDS(cd8_cite, "data/cd8_cite_analyzed.rds")

cat("\n=====  =====\n")
cat("CD8⁺ T ：", ncol(cd8_cite), "\n")
cat("：", length(unique(cd8_cite$seurat_clusters)), "\n")
cat("\nHLA-DR/CD27 ：\n")
print(table(cd8_cite$HLA_DR_CD27_group))

cat("\n✅ ！\n")


plan("sequential")