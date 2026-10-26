library(ggplot2)
library(ggpubr)
library(patchwork)

morandi_colors <- list(
  group = c("E-DP" = "#8C5A5A",
            "E-Trp" = "#6B7B6B",
            "H-Trp" = "#7A8B99"),
  cloneType = c("Single" = "#D4CFC5",
                "Small" = "#B5C1B4",
                "Medium" = "#8FA3A6",
                "Large" = "#6B7B8F",
                "Hyperexpanded" = "#8C5A5A"),
  gradient = c(low = "#D4CFC5", mid = "#B5C1B4", high = "#8C5A5A")
)

theme_pub <- theme_minimal(base_size = 7, base_family = "Arial") +
  theme(

    plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(b = 3)),

    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = unit(0.05, "cm"),

    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",

    panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.3, fill = NA, color = "black"),

    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_rect(linewidth = 0.3, fill = "grey95", color = "black")
  )

theme_set(theme_pub)

score_features <- c("GBD_Pathogenic1", "CD8_Treg1", "Exhaustion_Full1",
                    "Stemness_Score1", "Costim_Score1", "Migration_Score1",
                    "TissueResident_Score1", "Glycolysis_Score1", "OXPHOS_Score1")

score_features <- intersect(score_features, colnames(seurat_filtered@meta.data))

score_labels <- c(
  "GBD_Pathogenic1" = "GBD Pathogenic",
  "CD8_Treg1" = "CD8+ Treg",
  "Exhaustion_Full1" = "Exhaustion",
  "Stemness_Score1" = "Stemness",
  "Costim_Score1" = "Costimulation",
  "Migration_Score1" = "Migration",
  "TissueResident_Score1" = "Tissue Resident",
  "Glycolysis_Score1" = "Glycolysis",
  "OXPHOS_Score1" = "OXPHOS"
)

p_scores <- VlnPlot(seurat_filtered,
                    features = score_features,
                    group.by = "group",
                    pt.size = 0,
                    cols = morandi_colors$group,
                    ncol = 3) &
  theme_pub &
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p_scores <- p_scores +
  plot_annotation(title = "Functional Signatures Across Groups",
                  theme = theme_pub)

ggsave("Fig_FunctionalScores_by_Group.pdf", plot = p_scores,
       width = 7.5, height = 6, device = cairo_pdf)

if (!require("ComplexHeatmap")) {
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(-2, 0, 2), c("#D4CFC5", "#B5C1B4", "#8C5A5A"))

ht <- Heatmap(avg_mat_scaled,
              name = "Expression",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 7, fontfamily = "sans"),
              column_names_gp = gpar(fontsize = 7, fontfamily = "sans"),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 7, fontfamily = "sans"),
                labels_gp = gpar(fontsize = 6, fontfamily = "sans")
              ),
              column_title = "Key Gene Expression: Group  Clone Type",
              column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
              border = "white",
              width = unit(6, "cm"),
              height = unit(8, "cm"))

pdf("Fig_Heatmap_KeyGenes_ComplexHeatmap.pdf", width = 8, height = 6)
draw(ht)
dev.off()

p_umap_group <- DimPlot(seurat_filtered,
                        group.by = "group",
                        cols = morandi_colors$group,
                        pt.size = 0.5) +
  labs(title = "Group Distribution") +
  theme_pub +
  theme(legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank())

p_umap_clone <- DimPlot(seurat_filtered,
                        group.by = "cloneType_real",
                        cols = morandi_colors$cloneType,
                        pt.size = 0.5) +
  labs(title = "Clonal Expansion") +
  theme_pub +
  theme(legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank())

p_umap_combined <- p_umap_group + p_umap_clone +
  plot_annotation(tag_levels = 'A') &
  theme_pub

ggsave("Fig_UMAP_Combined.pdf", plot = p_umap_combined,
       width = 7.5, height = 3.5, device = cairo_pdf)

deg_summary <- data.frame(
  Gene = c("FOXP3", "CTLA4", "TIGIT", "TCF7", "CD28", "HLA-DRA", "CD27", "CD8A"),
  Expansion = c(8.46, 2.60, 0.27, -3.49, -2.09, -1.12, -0.81, -2.95),
  Sorting = c(-2.99, 0.42, NA, 1.52, 2.78, 2.79, 1.74, 1.57)
) %>%
  pivot_longer(cols = c(Expansion, Sorting),
               names_to = "Comparison", values_to = "log2FC") %>%
  mutate(Comparison = factor(Comparison,
                             levels = c("Expansion", "Sorting"),
                             labels = c("Expansion\n(E-Trp vs H-Trp)",
                                        "Sorting Timing\n(E-DP vs E-Trp)")))

p_bubble <- ggplot(deg_summary, aes(x = Comparison, y = Gene,
                                    size = abs(log2FC), fill = log2FC)) +
  geom_point(shape = 21, color = "white", stroke = 0.3) +
  scale_fill_gradient2(low = "#7A8B99", mid = "#D4CFC5", high = "#8C5A5A",
                       midpoint = 0, limits = c(-4, 9), na.value = "grey80") +
  scale_size_continuous(range = c(3, 12), guide = "none") +
  labs(x = "", y = "", fill = "log2FC") +
  ggtitle("Functional Remodeling of CD8+ T Cells") +
  theme_pub +
  theme(axis.text.x = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 7, face = "italic"),
        legend.key.height = unit(0.8, "cm"))

ggsave("Fig_Bubble_Remodeling.pdf", plot = p_bubble,
       width = 4.5, height = 4, device = cairo_pdf)

clone_summary <- seurat_filtered@meta.data %>%
  filter(!is.na(cloneType_real)) %>%
  group_by(group, cloneType_real) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count) * 100)

p_bar <- ggplot(clone_summary, aes(x = group, y = proportion, fill = cloneType_real)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = morandi_colors$cloneType) +
  labs(x = "", y = "Proportion (%)", fill = "Clone Type") +
  ggtitle("Clonal Composition by Group") +
  theme_pub

ggsave("Fig_ClonalProportion_Bar.pdf", plot = p_bar,
       width = 4, height = 3.5, device = cairo_pdf)

layout <- "
AABB
CCDD
EEFF
"

p_final <- wrap_plots(
  A = p_umap_group,
  B = p_umap_clone,
  C = p_bar,
  D = p_bubble,
  E = p_scores[[1]],
  F = p_scores[[2]],
  design = layout
) + plot_annotation(tag_levels = 'A') &
  theme_pub

ggsave("Fig_Comprehensive_Main.pdf", plot = p_final,
       width = 8.5, height = 10, device = cairo_pdf)

library(scRepertoire)
library(pheatmap)

overlap_matrix <- clonalOverlap(combined_tcr,
                                cloneCall = "gene",
                                method = "overlap")

heat_colors <- colorRampPalette(c("#D4CFC5", "#B5C1B4", "#6B7B8F", "#8C5A5A"))(100)

pheatmap(overlap_matrix,
         color = heat_colors,
         display_numbers = TRUE,
         number_format = "%.0f",
         number_color = "white",
         fontsize = 8,
         fontfamily = "sans",
         border_color = "white",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Clonal Overlap Between Samples",
         filename = "Fig_ClonalOverlap_Heatmap.pdf",
         width = 6, height = 5)

group_clones <- list(
  `E-DP` = unique(unlist(lapply(combined_tcr[grep("E-CD8DP", names(combined_tcr))],
                                function(x) x$CTgene))),
  `E-Trp` = unique(unlist(lapply(combined_tcr[grep("E-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene))),
  `H-Trp` = unique(unlist(lapply(combined_tcr[grep("H-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene)))
)

group_overlap <- matrix(0, nrow = 3, ncol = 3)
rownames(group_overlap) <- names(group_clones)
colnames(group_overlap) <- names(group_clones)

for (i in 1:3) {
  for (j in 1:3) {
    group_overlap[i, j] <- length(intersect(group_clones[[i]], group_clones[[j]]))
  }
}

cat("\n==========  ==========\n")
print(group_overlap)

jaccard_matrix <- matrix(0, nrow = 3, ncol = 3)
rownames(jaccard_matrix) <- names(group_clones)
colnames(jaccard_matrix) <- names(group_clones)

for (i in 1:3) {
  for (j in 1:3) {
    intersection <- length(intersect(group_clones[[i]], group_clones[[j]]))
    union_size <- length(union(group_clones[[i]], group_clones[[j]]))
    jaccard_matrix[i, j] <- intersection / union_size
  }
}

cat("\n========== Jaccard  ==========\n")
print(round(jaccard_matrix, 3))

pheatmap(jaccard_matrix,
         color = heat_colors,
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "white",
         fontsize = 10,
         fontfamily = "sans",
         border_color = "white",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Clonal Similarity (Jaccard Index)",
         filename = "Fig_JaccardSimilarity_Heatmap.pdf",
         width = 5, height = 4)

group_clones <- list(
  `E-DP` = unique(unlist(lapply(combined_tcr[grep("E-CD8DP", names(combined_tcr))],
                                function(x) x$CTgene))),
  `E-Trp` = unique(unlist(lapply(combined_tcr[grep("E-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene))),
  `H-Trp` = unique(unlist(lapply(combined_tcr[grep("H-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene)))
)

jaccard_matrix <- matrix(0, nrow = 3, ncol = 3,
                         dimnames = list(names(group_clones), names(group_clones)))

for (i in 1:3) {
  for (j in 1:3) {
    inter <- length(intersect(group_clones[[i]], group_clones[[j]]))
    uni <- length(union(group_clones[[i]], group_clones[[j]]))
    jaccard_matrix[i, j] <- inter / uni
  }
}

cat("\n========== Jaccard Similarity Between Groups ==========\n")
print(round(jaccard_matrix, 3))

jaccard_melt <- melt(jaccard_matrix)
colnames(jaccard_melt) <- c("Group1", "Group2", "Jaccard")

ggplot(jaccard_melt, aes(x = Group1, y = Group2, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(Jaccard, 3)), size = 3, color = "white", fontface = "bold") +
  scale_fill_gradientn(colors = c("#D4CFC5", "#B5C1B4", "#6B7B8F", "#8C5A5A")) +
  coord_fixed() +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "", fill = "Jaccard\nIndex") +
  ggtitle("Clonal Similarity (Jaccard) by Group")

ggsave("Fig_JaccardSimilarity_Group.pdf", width = 5, height = 4, device = cairo_pdf)

group_clones <- list(
  `E-DP` = unique(unlist(lapply(combined_tcr[grep("E-CD8DP", names(combined_tcr))],
                                function(x) x$CTgene))),
  `E-Trp` = unique(unlist(lapply(combined_tcr[grep("E-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene))),
  `H-Trp` = unique(unlist(lapply(combined_tcr[grep("H-CD8Trp", names(combined_tcr))],
                                 function(x) x$CTgene)))
)

jaccard_matrix <- matrix(0, nrow = 3, ncol = 3,
                         dimnames = list(names(group_clones), names(group_clones)))

for (i in 1:3) {
  for (j in 1:3) {
    inter <- length(intersect(group_clones[[i]], group_clones[[j]]))
    uni <- length(union(group_clones[[i]], group_clones[[j]]))
    jaccard_matrix[i, j] <- inter / uni
  }
}

cat("\n========== Jaccard Similarity Between Groups ==========\n")
print(round(jaccard_matrix, 3))

jaccard_melt <- melt(jaccard_matrix)
colnames(jaccard_melt) <- c("Group1", "Group2", "Jaccard")

ggplot(jaccard_melt, aes(x = Group1, y = Group2, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(Jaccard, 3)), size = 3, color = "white", fontface = "bold") +
  scale_fill_gradientn(colors = c("#D4CFC5", "#B5C1B4", "#6B7B8F", "#8C5A5A")) +
  coord_fixed() +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "", fill = "Jaccard\nIndex") +
  ggtitle("Clonal Similarity (Jaccard) by Group")

ggsave("Fig_JaccardSimilarity_Group.pdf", width = 5, height = 4, device = cairo_pdf)

p_rarefaction <- clonalRarefaction(combined_tcr,
                                   cloneCall = "gene",
                                   group.by = "sample") +
  scale_color_manual(values = c(
    "E-CD8DP-1" = "#8C5A5A", "E-CD8DP-2" = "#A87575",
    "E-CD8Trp-1" = "#6B7B6B", "E-CD8Trp-2" = "#8A9A8A",
    "H-CD8Trp-1" = "#7A8B99", "H-CD8Trp-2" = "#9AABBA"
  )) +
  theme_pub +
  labs(x = "Number of Cells Sampled", y = "Unique Clonotypes") +
  ggtitle("Clonal Rarefaction Curves")

ggsave("Fig_ClonalRarefaction.pdf", plot = p_rarefaction,
       width = 6, height = 4.5, device = cairo_pdf)

group_combined_tcr <- list(
  `E-DP` = do.call(rbind, combined_tcr[grep("E-CD8DP", names(combined_tcr))]),
  `E-Trp` = do.call(rbind, combined_tcr[grep("E-CD8Trp", names(combined_tcr))]),
  `H-Trp` = do.call(rbind, combined_tcr[grep("H-CD8Trp", names(combined_tcr))])
)

diversity_indices <- data.frame(
  Group = names(group_clones),
  Unique_Clones = sapply(group_clones, length),
  Total_Cells = sapply(group_combined_tcr, nrow)
) %>%
  mutate(
    Shannon = sapply(group_combined_tcr, function(x) {
      freq <- table(x$CTgene)
      -sum((freq/sum(freq)) * log(freq/sum(freq)))
    }),
    Simpson = sapply(group_combined_tcr, function(x) {
      freq <- table(x$CTgene)
      sum((freq/sum(freq))^2)
    }),
    Inverse_Simpson = 1 / Simpson,
    Chao1 = sapply(group_clones, function(x) {

      n_singletons <- sum(table(x) == 1)
      n_doubletons <- sum(table(x) == 2)
      length(x) + (n_singletons^2) / (2 * n_doubletons + 1)
    })
  )

cat("\n==========  ==========\n")
print(diversity_indices)

write.csv(diversity_indices, "Table_ClonalDiversity_Indices.csv", row.names = FALSE)

library(tidyr)
div_long <- diversity_indices %>%
  pivot_longer(cols = c(Unique_Clones, Shannon, Inverse_Simpson, Chao1),
               names_to = "Index", values_to = "Value")

ggplot(div_long, aes(x = Group, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = morandi_colors$group) +
  facet_wrap(~Index, scales = "free_y", nrow = 1) +
  labs(x = "", y = "Value") +
  ggtitle("Clonal Diversity Metrics by Group") +
  theme_pub +
  theme(legend.position = "none")

ggsave("Fig_DiversityMetrics_Bar.pdf", width = 8, height = 3.5, device = cairo_pdf)

v_usage_list <- list()
for (s in names(combined_tcr)) {
  df <- combined_tcr[[s]]

  v_alpha <- gsub("\\..*", "", df$TCR1)
  v_beta  <- gsub("\\..*", "", df$TCR2)
  v_all <- c(v_alpha, v_beta)
  v_all <- v_all[!is.na(v_all) & v_all != ""]
  v_usage_list[[s]] <- table(v_all)
}

v_usage_df <- do.call(rbind, lapply(names(v_usage_list), function(s) {
  data.frame(sample = s, V_gene = names(v_usage_list[[s]]),
             count = as.numeric(v_usage_list[[s]]), stringsAsFactors = FALSE)
})) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100)

p_v_usage_manual <- ggplot(v_usage_df, aes(x = V_gene, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(
    "E-CD8DP-1" = "#8C5A5A", "E-CD8DP-2" = "#A87575",
    "E-CD8Trp-1" = "#6B7B6B", "E-CD8Trp-2" = "#8A9A8A",
    "H-CD8Trp-1" = "#7A8B99", "H-CD8Trp-2" = "#9AABBA"
  )) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "right") +
  labs(x = "V Gene", y = "Proportion (%)") +
  ggtitle("TCR V Gene Usage by Sample")

ggsave("Fig_VGene_Usage_manual.pdf", plot = p_v_usage_manual, width = 12, height = 5, device = cairo_pdf)

j_usage_list <- list()
for (s in names(combined_tcr)) {
  df <- combined_tcr[[s]]

  j_alpha <- gsub(".*\\.", "", gsub("\\..*", "", df$TCR1))
  j_beta  <- gsub(".*\\.", "", gsub("\\..*", "", df$TCR2))
  j_all <- c(j_alpha, j_beta)
  j_all <- j_all[!is.na(j_all) & j_all != ""]
  j_usage_list[[s]] <- table(j_all)
}

j_usage_df <- do.call(rbind, lapply(names(j_usage_list), function(s) {
  data.frame(sample = s, J_gene = names(j_usage_list[[s]]),
             count = as.numeric(j_usage_list[[s]]), stringsAsFactors = FALSE)
})) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100)

p_j_usage_manual <- ggplot(j_usage_df, aes(x = J_gene, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(
    "E-CD8DP-1" = "#8C5A5A", "E-CD8DP-2" = "#A87575",
    "E-CD8Trp-1" = "#6B7B6B", "E-CD8Trp-2" = "#8A9A8A",
    "H-CD8Trp-1" = "#7A8B99", "H-CD8Trp-2" = "#9AABBA"
  )) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(x = "J Gene", y = "Proportion (%)") +
  ggtitle("TCR J Gene Usage by Sample")

ggsave("Fig_JGene_Usage_manual.pdf", plot = p_j_usage_manual, width = 12, height = 5, device = cairo_pdf)

clone_v_data <- seurat_filtered@meta.data %>%
  filter(!is.na(CTgene) & !is.na(cloneType_real)) %>%
  select(CTgene, cloneType_real, cloneSize_real, group) %>%
  distinct() %>%
  left_join(
    do.call(rbind, lapply(names(combined_tcr), function(s) {
      df <- combined_tcr[[s]]
      data.frame(
        CTgene = df$CTgene,
        TRAV = gsub("\\..*", "", df$TCR1),
        TRBV = gsub("\\..*", "", df$TCR2),
        stringsAsFactors = FALSE
      )
    })) %>% distinct(CTgene, .keep_all = TRUE),
    by = "CTgene"
  )

hyperexpanded_v <- clone_v_data %>%
  filter(cloneType_real == "Hyperexpanded") %>%
  pivot_longer(cols = c(TRAV, TRBV), names_to = "Chain", values_to = "V_gene") %>%
  filter(!is.na(V_gene) & V_gene != "") %>%
  group_by(V_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100) %>%
  arrange(desc(count)) %>%
  head(15)

ggplot(hyperexpanded_v, aes(x = reorder(V_gene, proportion), y = proportion)) +
  geom_bar(stat = "identity", fill = "#8C5A5A", width = 0.7) +
  coord_flip() +
  labs(x = "V Gene", y = "Proportion of Hyperexpanded Clones (%)") +
  ggtitle("V Gene Preference in Hyperexpanded Clones") +
  theme_pub

ggsave("Fig_Hyperexpanded_VGene_Preference.pdf", width = 6, height = 5, device = cairo_pdf)

library(dplyr)

v_data_list <- list()
for (s in names(combined_tcr)) {
  df <- combined_tcr[[s]]
  v_alpha <- gsub("\\..*", "", df$TCR1)
  v_beta  <- gsub("\\..*", "", df$TCR2)
  v_all <- c(v_alpha, v_beta)
  v_all <- v_all[!is.na(v_all) & v_all != ""]
  v_data_list[[s]] <- data.frame(sample = s, V_gene = v_all, stringsAsFactors = FALSE)
}
v_all_df <- do.call(rbind, v_data_list)
v_all_df$group <- case_when(
  grepl("^E-CD8DP", v_all_df$sample) ~ "E-DP",
  grepl("^E-CD8Trp", v_all_df$sample) ~ "E-Trp",
  grepl("^H-CD8Trp", v_all_df$sample) ~ "H-Trp"
)

top_v <- names(sort(table(v_all_df$V_gene), decreasing = TRUE)[1:15])

enrich_results <- list()
for (v in top_v) {
  groups <- c("E-DP", "E-Trp", "H-Trp")
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      g1 <- groups[i]; g2 <- groups[j]
      cnt1 <- sum(v_all_df$V_gene == v & v_all_df$group == g1)
      total1 <- sum(v_all_df$group == g1)
      cnt2 <- sum(v_all_df$V_gene == v & v_all_df$group == g2)
      total2 <- sum(v_all_df$group == g2)
      mat <- matrix(c(cnt1, total1 - cnt1, cnt2, total2 - cnt2), nrow = 2)
      p_val <- fisher.test(mat)$p.value
      enrich_results[[length(enrich_results) + 1]] <- data.frame(
        V_gene = v, Comparison = paste(g1, "vs", g2),
        Prop1 = round(cnt1/total1*100,2), Prop2 = round(cnt2/total2*100,2),
        P_value = p_val, stringsAsFactors = FALSE
      )
    }
  }
}
enrich_df <- do.call(rbind, enrich_results)
enrich_df$FDR <- p.adjust(enrich_df$P_value, method = "BH")

sig <- enrich_df[enrich_df$FDR < 0.05, ]
if (nrow(sig) > 0) {
  cat("\n==========  V  (FDR < 0.05) ==========\n")
  print(sig)
  write.csv(sig, "Table_VGene_Enrichment_Significant.csv", row.names = FALSE)
} else {
  cat("\n FDR < 0.05  V \n")
}

write.csv(enrich_df, "Table_VGene_Enrichment_Full.csv", row.names = FALSE)
cat("\n  V \n")