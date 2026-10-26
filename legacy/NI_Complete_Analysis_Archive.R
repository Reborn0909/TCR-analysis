








library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(grid)

setwd("DATA_ROOT")
out_base <- "DATA_ROOT"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 14))


cat("\n========== 1.  ==========\n")
cd8 <- readRDS("data/cd8_cite_final_ref.rds")
cd4 <- readRDS("data/cd4_cells_final_ref.rds")

if(!"final_subset" %in% colnames(cd8@meta.data)) stop("cd8final_subset")
if(!"final_subset" %in% colnames(cd4@meta.data)) stop("cd4final_subset")

cat("CD8:", ncol(cd8), " CD4:", ncol(cd4), "\n")


cat("\n========== 2.  ==========\n")
cd8$time_group <- ifelse(as.character(cd8$timepoint) == "Day0", "Day0", "Post")
cd4$time_group <- ifelse(as.character(cd4$timepoint) == "Day0", "Day0", "Post")

cd8$final_subset_new <- as.character(cd8$final_subset)
cd8$final_subset_new[cd8$final_subset %in% c("TEM","Teff (Temra)")] <- "TEM"
cd4$final_subset_new <- as.character(cd4$final_subset)
cd4$final_subset_new[cd4$final_subset == "Teff (Th1-like)*"] <- "Teff"

cat("CD8:", paste(unique(cd8$final_subset_new), collapse=", "), "\n")
cat("CD4:", paste(unique(cd4$final_subset_new), collapse=", "), "\n")


cat("\n========== 3.  ==========\n")
process_obj <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  if(!"percent.mt" %in% colnames(obj@meta.data))
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)
}

cd8 <- process_obj(cd8)
cd4 <- process_obj(cd4)
cat("\n")


cat("\n========== 4.  ==========\n")


cd8_markers <- c("TCF7","LEF1","SELL","CCR7","CD27","CD28",
                 "GZMH","GZMB","GNLY","NKG7","PRF1","FGFBP2",
                 "KLRG1","CX3CR1","EOMES","TBX21",
                 "PDCD1","TIGIT","HAVCR2","LAG3")

cd4_markers <- c("TCF7","LEF1","SELL","CCR7","CD27","CD28",
                 "GZMH","GZMB","GZMA","GZMK","PRF1",
                 "IFNG","TBX21","EOMES","KLRG1","CX3CR1",
                 "FOXP3","CTLA4","IL2RA","IKZF2",
                 "PDCD1","TIGIT","HAVCR2","LAG3","TOX")

hla_cp <- c("HLA-DRA","HLA-DRB1","CD27","CD28","ICOS",
            "PDCD1","TIGIT","CTLA4","IKZF2","FOXP3","TCF7")

reg_genes <- c("IL10","TGFB1","CTLA4","PDCD1","TIGIT","LAG3",
               "HAVCR2","ENTPD1","NT5E","CD274","FOXP3","PRDM1","TOX")


cols_cd8 <- c("Naive (Tn)"="#2E75B6","TCM"="#D95F02","TEM"="#7570B3")
cols_cd4 <- c("Naive (Tn)"="#2E75B6","TCM"="#D95F02","TEM"="#7570B3",
              "Teff"="#1B9E77","Treg"="#E7298A")
qc_colors <- c("Day0"="#2E75B6","Post"="#D95F02")
fp_low <- "lightgrey"
fp_high_blue <- "#2E75B6"
fp_high_red <- "#D95F02"
h_low <- "#2166AC"
h_mid <- "#F7F7F7"
h_high <- "#B2182B"
dp_colors <- c("HLA-DR+CD27+"="red","Other"="grey80")


cat("\n========== 5.  ==========\n")

plot_heatmap_ggplot <- function(mat, title = "") {
  mat <- t(scale(t(mat)))
  df <- as.data.frame(mat) %>% rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Group", values_to = "Expression")
  df$Gene <- factor(df$Gene, levels = rev(rownames(mat)))
  df$Group <- factor(df$Group, levels = colnames(mat))
  ggplot(df, aes(Group, Gene, fill = Expression)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = h_low, mid = h_mid, high = h_high, midpoint = 0) +
    labs(title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 11),
          plot.title = element_text(hjust = 0.5, size = 16),
          panel.grid = element_blank())
}

plot_combined_diff_heatmap <- function(fc_list, cell_type, title_suffix = "Post vs Day0") {
  all_genes <- unique(unlist(lapply(fc_list, names)))
  mat <- matrix(0, nrow = length(all_genes), ncol = length(fc_list))
  rownames(mat) <- all_genes
  colnames(mat) <- names(fc_list)
  for(i in seq_along(fc_list)) {
    mat[names(fc_list[[i]]), i] <- fc_list[[i]]
  }
  df <- as.data.frame(mat) %>% rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Comparison", values_to = "log2FC")
  df$Gene <- factor(df$Gene, levels = rev(rownames(mat)))
  df$Comparison <- factor(df$Comparison, levels = colnames(mat))
  ggplot(df, aes(Comparison, Gene, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", log2FC)), size = 3.5) +
    scale_fill_gradient2(low = h_low, mid = h_mid, high = h_high, midpoint = 0) +
    labs(title = paste(cell_type, ":", title_suffix)) +
    theme(axis.text.x = element_text(size = 13, face = "bold"),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14),
          panel.grid = element_blank())
}


cat("\n========== 6.  ==========\n")

report_ni <- function(obj, cell_type, marker_genes, cols_sub) {
  DefaultAssay(obj) <- "RNA"
  genes_use <- intersect(marker_genes, rownames(obj))


  cat(sprintf("  - %s ...\n", cell_type))
  print(VlnPlot(obj, "nFeature_RNA", group.by = "time_group", pt.size = 0, cols = qc_colors) +
          labs(title = paste(cell_type, "QC: nFeature_RNA")) + theme(legend.position = "none"))
  print(VlnPlot(obj, "nCount_RNA", group.by = "time_group", pt.size = 0, cols = qc_colors) +
          labs(title = paste(cell_type, "QC: nCount_RNA")) + theme(legend.position = "none"))
  print(VlnPlot(obj, "percent.mt", group.by = "time_group", pt.size = 0, cols = qc_colors) +
          labs(title = paste(cell_type, "QC: % MT")) + theme(legend.position = "none"))
  print(obj@meta.data %>% count(time_group) %>%
          ggplot(aes(time_group, n, fill = time_group)) + geom_bar(stat = "identity") +
          scale_fill_manual(values = qc_colors) + labs(title = paste(cell_type, ": Cell Counts")) +
          theme(legend.position = "none"))


  cat(sprintf("  - %s ...\n", cell_type))
  print(DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
          labs(title = paste(cell_type, ": Original Clusters (All)")))
  print(DimPlot(obj, group.by = "seurat_clusters", label = TRUE, split.by = "time_group", pt.size = 0.3) +
          labs(title = paste(cell_type, ": Original Clusters (Day0 vs Post)")))
  print(DotPlot(obj, features = genes_use, group.by = "seurat_clusters", dot.scale = 8) +
          labs(title = paste(cell_type, ": Original Cluster DotPlot")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)))


  cat(sprintf("  - %s ...\n", cell_type))
  print(DimPlot(obj, group.by = "final_subset_new", cols = cols_sub, pt.size = 0.3) +
          labs(title = paste(cell_type, ": Final Annotation (All)")))
  print(DimPlot(obj, group.by = "final_subset_new", cols = cols_sub, split.by = "time_group", pt.size = 0.3) +
          labs(title = paste(cell_type, ": Final Annotation (Day0 vs Post)")))
  print(DotPlot(obj, features = genes_use, group.by = "final_subset_new", dot.scale = 10) +
          labs(title = paste(cell_type, ": Final Annotation DotPlot")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)))


  cat(sprintf("  - %s ...\n", cell_type))
  prop <- obj@meta.data %>% count(time_group, final_subset_new) %>%
    group_by(time_group) %>% mutate(pct = n/sum(n)*100)
  print(ggplot(prop, aes(time_group, pct, fill = final_subset_new)) +
          geom_bar(stat = "identity", position = "fill", width = 0.7) +
          scale_fill_manual(values = cols_sub) +
          labs(title = paste(cell_type, ": Subset Proportion"), y = "Proportion") +
          theme(axis.text.x = element_text(size = 14)))


  cat(sprintf("  - %s ...\n", cell_type))
  key_blue <- intersect(c("TCF7","SELL","GZMH","GZMB","KLRG1","PDCD1","TIGIT","FOXP3"), rownames(obj))
  if(length(key_blue) > 0) {
    pb <- lapply(key_blue, function(g) {
      FeaturePlot(obj, g, order = TRUE, cols = c(fp_low, fp_high_blue), pt.size = 0.2) +
        labs(title = g) + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10))
    })
    print(wrap_plots(pb, ncol = 4))
  }

  cp_genes <- intersect(hla_cp, rownames(obj))
  if(length(cp_genes) > 0) {
    pr <- lapply(cp_genes, function(g) {
      FeaturePlot(obj, g, order = TRUE, cols = c(fp_low, fp_high_red), pt.size = 0.2) +
        labs(title = g) + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10))
    })
    print(wrap_plots(pr, ncol = 5))
  }


  cat(sprintf("  - %s ...\n", cell_type))
  vln_groups <- list(
    "Trafficking/Homing"      = c("CCR7","SELL","CD28","IL7R"),
    "Stemness/Memory"         = c("TCF7","LEF1","CD27"),
    "Effector/Cytotoxicity"   = c("GZMH","GZMB","GNLY","NKG7","PRF1","FGFBP2"),
    "Terminal Differentiation"= c("KLRG1","CX3CR1","EOMES","TBX21"),
    "Exhaustion/Senescence"   = c("PDCD1","TIGIT","HAVCR2","LAG3")
  )
  if(cell_type == "CD4") {
    vln_groups[["Treg Core"]] <- c("FOXP3","CTLA4","IL2RA","IKZF2")
  }
  for(grp in names(vln_groups)) {
    genes_use_vln <- intersect(vln_groups[[grp]], rownames(obj))
    if(length(genes_use_vln) > 0) {
      print(VlnPlot(obj, features = genes_use_vln, group.by = "final_subset_new",
                    cols = cols_sub, ncol = length(genes_use_vln), pt.size = 0) +
              plot_annotation(title = paste(cell_type, ":", grp),
                              theme = theme(plot.title = element_text(size = 14))))
    }
  }


  cat(sprintf("  - %s ...\n", cell_type))
  if(length(genes_use) > 0) {
    avg <- AverageExpression(obj, features = genes_use,
                             group.by = "final_subset_new", assays = "RNA")[["RNA"]]
    if(ncol(avg) > 1 && nrow(avg) > 1) {
      print(plot_heatmap_ggplot(avg, paste(cell_type, ": Average Expression")))
    }
  }


  cat(sprintf("  - %s ...\n", cell_type))
  fc_list <- list()
  for(sub in names(cols_sub)) {
    if(!sub %in% unique(obj$final_subset_new)) next
    obj_sub <- obj[, obj$final_subset_new == sub]
    if(length(unique(obj_sub$time_group)) < 2) next
    avg_tp <- tryCatch({
      AverageExpression(obj_sub, features = genes_use,
                        group.by = "time_group", assays = "RNA")[["RNA"]]
    }, error = function(e) NULL)
    if(is.null(avg_tp)) next
    if("Day0" %in% colnames(avg_tp) & "Post" %in% colnames(avg_tp) & nrow(avg_tp) > 0) {
      fc_vec <- log2((avg_tp[,"Post"] + 0.01) / (avg_tp[,"Day0"] + 0.01))
      fc_list[[sub]] <- fc_vec
    }
  }
  if(length(fc_list) >= 2) {
    print(plot_combined_diff_heatmap(fc_list, cell_type, "Post vs Day0 per Subset"))
  }
}


pdf_base <- file.path(out_base, "NI_Base_Report.pdf")
cat(":", pdf_base, "\n")
pdf(pdf_base, width = 18, height = 14, onefile = TRUE)
report_ni(cd8, "CD8", cd8_markers, cols_cd8)
report_ni(cd4, "CD4", cd4_markers, cols_cd4)
dev.off()
cat("✅ \n")


cat("\n========== 7. DP ==========\n")


calculate_q25_threshold <- function(obj, gene) {
  expr <- FetchData(obj, gene)[,1]
  quantile(expr[expr > 0], 0.25)
}

th_hla_cd8 <- calculate_q25_threshold(cd8, "HLA-DRA")
th_cd27_cd8 <- calculate_q25_threshold(cd8, "CD27")
th_hla_cd4 <- calculate_q25_threshold(cd4, "HLA-DRA")
th_cd27_cd4 <- calculate_q25_threshold(cd4, "CD27")

cat(sprintf("CD8: HLA-DRA=%.3f CD27=%.3f\n", th_hla_cd8, th_cd27_cd8))
cat(sprintf("CD4: HLA-DRA=%.3f CD27=%.3f\n", th_hla_cd4, th_cd27_cd4))

cd8$DP <- ifelse(FetchData(cd8,"HLA-DRA")[,1] > th_hla_cd8 &
                   FetchData(cd8,"CD27")[,1] > th_cd27_cd8, "HLA-DR+CD27+", "Other")
cd4$DP <- ifelse(FetchData(cd4,"HLA-DRA")[,1] > th_hla_cd4 &
                   FetchData(cd4,"CD27")[,1] > th_cd27_cd4, "HLA-DR+CD27+", "Other")

cat("CD8 DP:", round(mean(cd8$DP=="HLA-DR+CD27+")*100, 1), "%\n")
cat("CD4 DP:", round(mean(cd4$DP=="HLA-DR+CD27+")*100, 1), "%\n")


reg_genes_cd8 <- intersect(reg_genes, rownames(cd8))
reg_genes_cd4 <- intersect(reg_genes, rownames(cd4))
cd8 <- AddModuleScore(cd8, list(reg_genes_cd8), name="RegScore")
cd4 <- AddModuleScore(cd4, list(reg_genes_cd4), name="RegScore")


cd8_nn <- subset(cd8, final_subset_new!="Naive (Tn)")
cd8_nn$Group <- factor(ifelse(cd8_nn$DP=="HLA-DR+CD27+","DP",
                              as.character(cd8_nn$final_subset_new)),
                       levels=c("DP","TCM","TEM"))
cd4_nn <- subset(cd4, final_subset_new!="Naive (Tn)")
cd4_nn$Group <- factor(ifelse(cd4_nn$DP=="HLA-DR+CD27+","DP",
                              as.character(cd4_nn$final_subset_new)),
                       levels=c("DP","TCM","TEM","Teff","Treg"))


wt1_cd8 <- wilcox.test(RegScore1~Group,
                       data=cd8_nn@meta.data[cd8_nn$Group%in%c("DP","TCM"),],
                       alternative="greater")
wt1_cd4 <- wilcox.test(RegScore1~Group,
                       data=cd4_nn@meta.data[cd4_nn$Group%in%c("DP","TCM"),],
                       alternative="greater")
cd8_tcm <- subset(cd8, final_subset_new=="TCM")
cd4_tcm <- subset(cd4, final_subset_new=="TCM")
wt2_cd8 <- wilcox.test(RegScore1~DP, data=cd8_tcm@meta.data, alternative="greater")
wt2_cd4 <- wilcox.test(RegScore1~DP, data=cd4_tcm@meta.data, alternative="greater")


pdf_dp <- file.path(out_base, "NI_DP_Report.pdf")
cat("DP:", pdf_dp, "\n")
pdf(pdf_dp, width = 18, height = 14)


print(wrap_plots(
  DimPlot(cd8, group.by="DP", cols=dp_colors, pt.size=0.4) + ggtitle("CD8: DP Cells") + NoLegend(),
  DimPlot(cd4, group.by="DP", cols=dp_colors, pt.size=0.4) + ggtitle("CD4: DP Cells") + NoLegend(),
  ncol=2))


print(VlnPlot(cd8_nn, "RegScore1", group.by="Group", pt.size=0,
              cols=c("DP"="red","TCM"="#D95F02","TEM"="#7570B3")) +
        ggtitle("CD8 non-Naive: Regulatory Score"))
print(VlnPlot(cd4_nn, "RegScore1", group.by="Group", pt.size=0,
              cols=c("DP"="red","TCM"="#D95F02","TEM"="#7570B3","Teff"="#1B9E77","Treg"="#E7298A")) +
        ggtitle("CD4 non-Naive: Regulatory Score"))


if(length(reg_genes_cd8) > 0) {
  avg_cd8 <- AverageExpression(cd8_nn, features = reg_genes_cd8,
                               group.by = "Group", assays = "RNA")[["RNA"]]
  if(ncol(avg_cd8) > 1 && nrow(avg_cd8) > 1) {
    pheatmap(t(scale(t(avg_cd8))),
             main = "CD8: Regulatory Genes (non-Naive)",
             cluster_cols = FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
  }
}
if(length(reg_genes_cd4) > 0) {
  avg_cd4 <- AverageExpression(cd4_nn, features = reg_genes_cd4,
                               group.by = "Group", assays = "RNA")[["RNA"]]
  if(ncol(avg_cd4) > 1 && nrow(avg_cd4) > 1) {
    pheatmap(t(scale(t(avg_cd4))),
             main = "CD4: Regulatory Genes (non-Naive)",
             cluster_cols = FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
  }
}


print(wrap_plots(
  VlnPlot(cd8_tcm, "RegScore1", group.by="DP", pt.size=0, cols=dp_colors) +
    ggtitle("CD8 TCM: DP vs non-DP"),
  VlnPlot(cd4_tcm, "RegScore1", group.by="DP", pt.size=0, cols=dp_colors) +
    ggtitle("CD4 TCM: DP vs non-DP"),
  ncol=2))


grid.newpage()
grid.text(
  sprintf("Wilcoxon Rank-Sum Test (one-sided)\n\n%s\n%s\n\n%s\n%s",
          paste0("CD8 DP vs TCM: p = ", format(wt1_cd8$p.value, digits=3, scientific=TRUE)),
          paste0("CD4 DP vs TCM: p = ", format(wt1_cd4$p.value, digits=3, scientific=TRUE)),
          paste0("CD8 TCM DP vs non-DP: p = ", format(wt2_cd8$p.value, digits=3, scientific=TRUE)),
          paste0("CD4 TCM DP vs non-DP: p = ", format(wt2_cd4$p.value, digits=3, scientific=TRUE))),
  x = 0.05, y = 0.5, just = "left",
  gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()
cat("✅ DP\n")


cat("\n========== 8.  ==========\n")


cd8$C5_DP <- ifelse(cd8$seurat_clusters==5, "C5 (DP)", as.character(cd8$final_subset_new))
cd8_nn2 <- subset(cd8, final_subset_new!="Naive (Tn)")
cd8_nn2$Group <- factor(ifelse(cd8_nn2$seurat_clusters==5, "C5 (DP)",
                               as.character(cd8_nn2$final_subset_new)),
                        levels=c("C5 (DP)","TCM","TEM"))
wt_c5 <- wilcox.test(RegScore1~Group,
                     data=cd8_nn2@meta.data[cd8_nn2$Group%in%c("C5 (DP)","TCM"),],
                     alternative="greater")


cd4_tcm_sub <- subset(cd4, final_subset_new=="TCM")
cd4_tcm_sub$DP_tcm <- ifelse((FetchData(cd4_tcm_sub,"HLA-DRA")[,1] > 0 |
                                FetchData(cd4_tcm_sub,"HLA-DRB1")[,1] > 0) &
                               FetchData(cd4_tcm_sub,"CD27")[,1] > 0, "DP", "non-DP")
wt_tcm_dp <- wilcox.test(RegScore1~DP_tcm, data=cd4_tcm_sub@meta.data, alternative="greater")


pdf_sup <- file.path(out_base, "NI_Supplementary_Report.pdf")
cat(":", pdf_sup, "\n")
pdf(pdf_sup, width = 18, height = 14)

print(DimPlot(cd8_nn2, group.by="Group", pt.size=0.4,
              cols=c("C5 (DP)"="red","TCM"="#D95F02","TEM"="#7570B3")) +
        ggtitle("CD8: C5 (DP) in non-Naive"))
print(VlnPlot(cd8_nn2, "RegScore1", group.by="Group", pt.size=0,
              cols=c("C5 (DP)"="red","TCM"="#D95F02","TEM"="#7570B3")) +
        ggtitle("CD8: Regulatory Score by Group"))

print(DimPlot(cd4_tcm_sub, group.by="DP_tcm", pt.size=0.4,
              cols=c("DP"="red","non-DP"="grey60")) +
        ggtitle("CD4 TCM: DP Distribution"))
print(VlnPlot(cd4_tcm_sub, "RegScore1", group.by="DP_tcm", pt.size=0,
              cols=c("DP"="red","non-DP"="grey60")) +
        ggtitle("CD4 TCM: Regulatory Score by DP"))

grid.newpage()
grid.text(
  sprintf("Supplementary Statistical Results\n\n%s\n%s",
          paste0("CD8 C5 (DP) vs TCM: p = ", format(wt_c5$p.value, digits=3, scientific=TRUE)),
          paste0("CD4 TCM DP vs non-DP: p = ", format(wt_tcm_dp$p.value, digits=3, scientific=TRUE))),
  x = 0.05, y = 0.5, just = "left",
  gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()
cat("✅ \n")


cat("\n========== 9.  ==========\n")
saveRDS(cd8, file.path(out_base, "NI_cd8_final.rds"))
saveRDS(cd4, file.path(out_base, "NI_cd4_final.rds"))


cat("\n========== 10.  ==========\n")
sink(file.path(out_base, "NI_Analysis_Summary.txt"))
cat("========================================\n")
cat("NI\n")
cat("========================================\n")
cat(":", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
cat("CD8:", ncol(cd8), "\n")
cat("CD4:", ncol(cd4), "\n\n")
cat("CD8:\n")
print(table(cd8$final_subset_new))
cat("\nCD4:\n")
print(table(cd4$final_subset_new))
cat("\nCD8 DP:", round(mean(cd8$DP=="HLA-DR+CD27+")*100, 2), "%\n")
cat("CD4 DP:", round(mean(cd4$DP=="HLA-DR+CD27+")*100, 2), "%\n")
cat("\n:\n")
cat("CD8 DP vs TCM p:", format(wt1_cd8$p.value, digits=3, scientific=TRUE), "\n")
cat("CD4 DP vs TCM p:", format(wt1_cd4$p.value, digits=3, scientific=TRUE), "\n")
cat("CD8 C5 vs TCM p:", format(wt_c5$p.value, digits=3, scientific=TRUE), "\n")
cat("\n:\n")
cat("1. NI_Base_Report.pdf - \n")
cat("2. NI_DP_Report.pdf - DP\n")
cat("3. NI_Supplementary_Report.pdf - \n")
cat("4. NI_cd8_final.rds - CD8\n")
cat("5. NI_cd4_final.rds - CD4\n")
cat("6. NI_Analysis_Summary.txt - \n")
cat("========================================\n")
sink()

cat("\n✅ NI！\n")
cat(":", out_base, "\n")
cat(":\n")
cat("  1. NI_Base_Report.pdf\n")
cat("  2. NI_DP_Report.pdf\n")
cat("  3. NI_Supplementary_Report.pdf\n")
cat("  4. NI_cd8_final.rds\n")
cat("  5. NI_cd4_final.rds\n")
cat("  6. NI_Analysis_Summary.txt\n")