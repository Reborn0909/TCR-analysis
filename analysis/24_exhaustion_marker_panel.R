




suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
setwd("DATA_ROOT")
OUT <- "figures/exhaustion"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


PANEL <- list(
  "Inhibitory receptor" = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","CD160","VSIR","CD244","LAIR1"),
  "Transcription factor" = c("TOX","TOX2","NR4A1","NR4A2","NR4A3","BATF","IRF4","PRDM1","ID2","EOMES","TBX21","IKZF2"),
  "Cytokine / chemokine" = c("CXCL13","CCL5","IFNG","GZMB","PRF1"),
  "Surface / other" = c("ENTPD1","ITGAE","TNFRSF9","TNFRSF18","SLAMF7","LAYN","CD27","CD28"),
  "Memory / progenitor (reference)" = c("TCF7","SELL","CCR7","IL7R","LEF1"))
ALLG <- unique(unlist(PANEL))

o <- readRDS("DATA_ROOT")
o$group <- factor(o$group, levels=c("H-Trp","E-Trp","E-DP"))
o$donor <- ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2")
cat("  :", ncol(o), "\n")

genes <- intersect(ALLG, rownames(o))
cat("  :", length(genes), "/", length(ALLG), "\n")
miss <- setdiff(ALLG, genes); if (length(miss)) cat("  :", paste(miss, collapse=", "), "\n")


emb <- FetchData(o, vars=c(genes, "group", "donor"), layer="data")
long <- emb %>% pivot_longer(all_of(genes), names_to="gene", values_to="expr") %>%
  group_by(donor, group, gene) %>%
  summarise(mean_log1p = mean(expr), pct_expressing = 100*mean(expr > 0), n_cells = n(), .groups="drop")


w <- long %>% select(donor, group, gene, mean_log1p) %>%
  pivot_wider(names_from=group, values_from=mean_log1p)
d <- w %>% mutate(delta_Pre_vs_SS = `E-Trp` - `H-Trp`,
                  delta_Post_vs_Pre = `E-DP` - `E-Trp`) %>%
  select(donor, gene, delta_Pre_vs_SS, delta_Post_vs_Pre)
agg <- d %>% group_by(gene) %>%
  summarise(d1 = delta_Pre_vs_SS[donor=="Donor 1"], d2 = delta_Pre_vs_SS[donor=="Donor 2"],
            mean_delta = mean(delta_Pre_vs_SS),
            up_both = all(delta_Pre_vs_SS > 0), down_both = all(delta_Pre_vs_SS < 0),
            .groups="drop") %>%
  mutate(category = sapply(gene, function(g) names(PANEL)[sapply(PANEL, function(x) g %in% x)][1])) %>%
  arrange(desc(mean_delta))

write.csv(agg, "results/exhaustion_marker_panel_donor_delta.csv", row.names=FALSE)
write.csv(long, "results/exhaustion_marker_panel_donor_expression.csv", row.names=FALSE)

cat("\n=== CD8 Trp (Pre-sel)  SS 、 ===\n")
up <- agg %>% filter(up_both) %>% arrange(desc(mean_delta))
for (i in seq_len(nrow(up))) cat(sprintf("  %-10s %-24s  D1 %+0.3f  D2 %+0.3f   %+0.3f\n",
  up$gene[i], up$category[i], up$d1[i], up$d2[i], up$mean_delta[i]))
cat("\n=== 、 ===\n")
dn <- agg %>% filter(down_both) %>% arrange(mean_delta)
for (i in seq_len(nrow(dn))) cat(sprintf("  %-10s %-24s  D1 %+0.3f  D2 %+0.3f   %+0.3f\n",
  dn$gene[i], dn$category[i], dn$d1[i], dn$d2[i], dn$mean_delta[i]))
cat("\n  ():", sum(agg$up_both), " :", sum(agg$down_both), " :", sum(!agg$up_both & !agg$down_both), "\n")
