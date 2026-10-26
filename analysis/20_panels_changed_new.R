





suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(grid)
  library(showtext); library(sysfonts)
})

font_add("msyh", regular = "C:/Windows/Fonts/msyh.ttc",
                  bold    = "C:/Windows/Fonts/msyhbd.ttc")
showtext_auto()
showtext_opts(dpi = 600)
setwd("DATA_ROOT")
RES <- "results"; OUT <- "figures/panels"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)


sample_colors <- c("E-CD8DP-1"="\u0023C8956B", "E-CD8DP-2"="\u0023D9AE8A",
                   "E-CD8Trp-1"="\u0023B05A5A", "E-CD8Trp-2"="\u0023C27A7A",
                   "H-CD8Trp-1"="\u0023A87575", "H-CD8Trp-2"="\u00238BA6B2")
sample_colors["H-CD8Trp-1"] <- "\u00236B8E9B"
group_colors  <- c("H-Trp"="\u00236B8E9B", "E-Trp"="\u0023B05A5A", "E-DP"="\u0023C8956B")


comp_colors   <- c("E-Trp vs H-Trp"="#6B8E9B", "E-DP vs E-Trp"="#B05A5A", "E-DP vs H-Trp"="#C8956B")
clone_colors  <- c("Single"="\u00236B8E9B", "Small"="\u00238BA6B2", "Medium"="\u0023C8956B",
                   "Large"="\u0023B05A5A", "Hyperexpanded"="\u00238C3A3A")
GRP <- c("H-Trp","E-Trp","E-DP")
SAMP <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")
DON  <- c("Donor1","Donor2")


theme_pub <- theme_minimal(base_size = 7, base_family = "msyh") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(b = 3)),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        axis.ticks.length = unit(0.05, "cm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"), legend.position = "right",
        panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 0.3, fill = NA, color = "black"),
        strip.text = element_text(size = 7, face = "bold"),
        strip.background = element_rect(linewidth = 0.3, fill = "grey95", color = "black"))
theme_set(theme_pub)

save_panel <- function(p, name, w = 3.4, h = 2.6) {
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h, family = "msyh"); print(p); dev.off()
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, type = "cairo")
  cat("  :", name, "\n")
}


qc <- read_csv(file.path(RES,"main_preQC_cell_counts.csv"), show_col_types = FALSE) %>%
  mutate(sample = factor(sample, levels = SAMP))
p1 <- qc %>% select(sample, pre_QC_cells, post_QC_cells) %>%
  pivot_longer(-sample, names_to="stage", values_to="cells") %>%
  mutate(stage = recode(stage, pre_QC_cells="QC ", post_QC_cells="QC ",
                        .default = stage),
         stage = factor(stage, levels=c("QC ","QC "))) %>%
  ggplot(aes(sample, cells, fill = sample, alpha = stage)) +
  geom_col(position = position_dodge(width=.75), width=.68, colour="grey30", linewidth=.2) +
  scale_fill_manual(values = sample_colors, guide="none") +
  scale_alpha_manual(values = c("QC "=.45, "QC "=1), name=NULL,
                     guide = guide_legend(override.aes = list(fill = "#8C5A5A"))) +
  scale_y_continuous(labels = comma, expand = expansion(mult=c(0,.08))) +
  labs(title = "P01  QC ", x = NULL, y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_panel(p1, "P01_QC_pre_post_cells", 4.2, 2.8)


attrb <- read_csv(file.path(RES,"main_raw_QC_threshold_attribution.csv"), show_col_types = FALSE)
names(attrb) <- c("sample","pre","post","feat_low","feat_high","umi_high","mt_high","hb_high")
THR <- c("percent.mt ≥ 10"="#C8956B", "nCount ≥ 80000"="#B05A5A",
         "nFeature ≤ 200"="#6B8E9B", "nFeature ≥ 10000"="#8BA6B2", "percent.HB ≥ 0.35"="#D4CFC5")
a2 <- attrb %>% mutate(sample = factor(sample, levels = SAMP)) %>%
  select(sample, `percent.mt ≥ 10`=mt_high, `nCount ≥ 80000`=umi_high,
         `nFeature ≤ 200`=feat_low, `nFeature ≥ 10000`=feat_high, `percent.HB ≥ 0.35`=hb_high) %>%
  pivot_longer(-sample, names_to="threshold", values_to="removed") %>%
  mutate(threshold = factor(threshold, levels = names(THR)))
p2 <- ggplot(a2, aes(sample, removed, fill = threshold)) +
  geom_col(width=.68, colour="grey30", linewidth=.2) +
  scale_fill_manual(values = THR, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult=c(0,.08))) +
  labs(title = "P02 ", x = NULL, y = "") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="right",
        legend.key.size = unit(.3,"cm"))
save_panel(p2, "P02_QC_threshold_attribution", 4.6, 2.8)


db <- read_csv("results/main_doublet_assessment.csv", show_col_types = FALSE) %>%
  mutate(sample = factor(sample, levels = SAMP))
p3 <- db %>% select(sample, `≥0.20`=frac_ge_0.20, `≥0.33`=frac_ge_0.33, `≥0.50`=frac_ge_0.50) %>%
  pivot_longer(-sample, names_to="cut", values_to="pct") %>%
  mutate(cut = factor(cut, levels=c("≥0.20","≥0.33","≥0.50"))) %>%
  ggplot(aes(sample, pct, fill = sample, alpha = cut)) +
  geom_col(position=position_dodge(width=.75), width=.66, colour="grey30", linewidth=.2) +
  geom_point(data = db, aes(sample, expected_rate_pct_10x, alpha=NULL), inherit.aes=FALSE,
             shape=4, size=1.5, colour="grey20", stroke=.6) +
  scale_fill_manual(values = sample_colors, guide="none") +
  scale_alpha_manual(values=c(.45,.75,1), name="",
                     guide = guide_legend(override.aes = list(fill = "#8C5A5A"))) +
  scale_y_continuous(expand = expansion(mult=c(0,.08))) +
  labs(title = "P03 （× = 10x ）", x = NULL, y = " (%)") +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top")
save_panel(p3, "P03_doublet_assessment", 4.4, 2.8)


ds <- read_csv(file.path(RES,"main_TCR_equal_cell_downsample_summary.csv"), show_col_types=FALSE) %>%
  filter(clonotype_definition == "CTstrict") %>%
  mutate(sample = factor(sample, levels = SAMP), donor = factor(donor, levels = DON),
         group = factor(group, levels = GRP))
p4 <- ggplot(ds, aes(group, downsample_median, colour = sample, group = sample)) +
  geom_errorbar(aes(ymin = sampling_lower_2.5pct, ymax = sampling_upper_97.5pct),
                width=.08, linewidth=.35) +
  geom_line(linewidth=.45) + geom_point(size=1.5) +
  facet_wrap(~ donor, nrow=1) +
  scale_colour_manual(values = sample_colors, name=NULL) +
  labs(title = "P04 TCR Shannon （n=6,752）", x=NULL, y=" Shannon ") +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
save_panel(p4, "P04_TCR_equal_cell_downsample", 5.2, 2.8)


dp <- read_csv(file.path(RES,"main_depth_matched_TCR_gene_level.csv"), show_col_types=FALSE) %>%
  mutate(sample = factor(sample, levels = SAMP), floor_umis = factor(floor_umis))
p5 <- ggplot(dp, aes(floor_umis, downsample_median, colour=sample, group=sample)) +
  geom_line(linewidth=.45) + geom_point(size=1.3) +
  scale_colour_manual(values = sample_colors, name=NULL) +
  labs(title = "P05 （VDJ UMI ）", x="UMI ≥", y=" Shannon ") +
  theme(legend.position="right")
save_panel(p5, "P05_depth_floor_sweep", 4.4, 2.8)


pb <- data.frame(
  comparison = factor(c("E-Trp vs H-Trp","E-DP vs E-Trp","E-DP vs H-Trp"),
                      levels=c("E-Trp vs H-Trp","E-DP vs E-Trp","E-DP vs H-Trp")),
  tested = 26083,
  fdr = c(12748, 15585, 18193),
  lfc = c(9808, 12241, 14980))
pb %>% select(comparison, `FDR<0.05`=fdr, ` |log2FC|>1`=lfc) %>%
  pivot_longer(-comparison, names_to="criterion", values_to="n") %>%
  ggplot(aes(comparison, n, fill=comparison, alpha=criterion)) +
  geom_col(position=position_dodge(width=.72), width=.64, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=26083, linewidth=.3, linetype="dashed", colour="grey40") +
  annotate("text", x=.62, y=26083, label=" 26,083", hjust=0, vjust=-.4, size=2) +
  scale_fill_manual(values=comp_colors, guide="none") +
  scale_alpha_manual(values=c(.5,1), name=NULL) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0,.12))) +
  labs(title = "P06  pseudobulk ", x=NULL, y="") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="top") -> p6
save_panel(p6, "P06_pseudobulk_significant", 4.6, 2.8)


sg <- read_csv(file.path(RES,"main_recomputed_signature_donor_summary.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP), donor=factor(donor, levels=DON),
         group=factor(group, levels=GRP))
p7 <- ggplot(sg, aes(group, mean_control_adjusted_score, fill=sample)) +
  geom_col(width=.66, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey30") +
  facet_grid(signature ~ donor, scales="free_y") +
  scale_fill_manual(values=sample_colors, name=NULL) +
  labs(title="P07 signature （， cell-level P）", x=NULL, y="") +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
save_panel(p7, "P07_signature_donor_level", 5.0, 3.0)


eg <- read_csv(file.path(RES,"main_requested_eight_genes_donor_expression.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP))
ord <- c("TCF7","CCR7","IL7R","LEF1","TOX","HAVCR2","ENTPD1","LAG3")
eg$gene <- factor(eg$gene, levels=rev(ord))
p8 <- ggplot(eg, aes(sample, gene, fill=mean_log1p_normalized_RNA)) +
  geom_tile(colour="white", linewidth=.4) +
  geom_text(aes(label=sprintf("%.2f", mean_log1p_normalized_RNA)), size=1.6, colour="grey15") +
  scale_fill_gradientn(colours=c("#D4CFC5", "#B5C1B4", "#8C5A5A"),
                       name="log1p") +
  labs(title="P08 memory / exhaustion （-）", x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())
save_panel(p8, "P08_eight_genes_heatmap", 4.4, 2.8)


sc <- read_csv(file.path(RES,"main_requested_memory_exhaustion_donor_scores.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP), donor=factor(donor, levels=DON),
         group=factor(group, levels=GRP))
p9 <- sc %>% select(sample, group, donor, memory4=memory4_control_adjusted_mean,
                    exhaustion4=exhaustion4_control_adjusted_mean) %>%
  pivot_longer(c(memory4, exhaustion4), names_to="score", values_to="value") %>%
  ggplot(aes(group, value, fill=sample)) +
  geom_col(width=.66, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey30") +
  facet_grid(score ~ donor, scales="free_y") +
  scale_fill_manual(values=sample_colors, name=NULL) +
  labs(title="P09 memory4  exhaustion4 （）", x=NULL, y="") +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
save_panel(p9, "P09_memory_exhaustion_scores", 5.0, 3.0)


tr <- read_csv(file.path(RES,"main_productive_TCR_donor_trends.csv"), show_col_types=FALSE) %>%
  mutate(donor=factor(donor, levels=DON))
p10 <- ggplot(tr, aes(comparison, downsample_delta, fill=donor)) +
  geom_col(position=position_dodge(width=.72), width=.64, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey30") +
  facet_wrap(~ definition, nrow=1, scales="free_y") +
  scale_fill_manual(values=c("Donor1"="#8C5A5A","Donor2"="#A87575"), name=NULL) +
  labs(title="P10 （8/8 ）", x=NULL, y=" Shannon ") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="top")
save_panel(p10, "P10_donor_direction", 5.0, 2.8)

cat("\n panel ", OUT, "\n")
cat(" PDF:", length(list.files(OUT, pattern="\\.pdf$")), "\n")
