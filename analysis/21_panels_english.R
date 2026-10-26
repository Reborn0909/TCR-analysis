




suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(scales); library(grid)
})
setwd(Sys.getenv("DATA_ROOT", "."))
RES <- "results"; OUT <- "figures/panels_en"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)


sample_colors <- c("E-CD8DP-1"="#C8956B", "E-CD8DP-2"="#D9AE8A",
                   "E-CD8Trp-1"="#B05A5A", "E-CD8Trp-2"="#C27A7A",
                   "H-CD8Trp-1"="#6B8E9B", "H-CD8Trp-2"="#8BA6B2")
group_colors  <- c("H-Trp"="#6B8E9B", "E-Trp"="#B05A5A", "E-DP"="#C8956B")
COND <- c("H-Trp"="SS", "E-Trp"="Pre-sel", "E-DP"="Post-sel")
GRP  <- c("H-Trp","E-Trp","E-DP")
SAMP <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")
DON  <- c("Donor 1","Donor 2")

theme_pub <- theme_minimal(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, margin = margin(b = 3)),
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

save_panel <- function(p, name, w = 3.6, h = 2.7) {
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, type = "cairo")
  cat("  saved:", name, "\n")
}


qc <- read_csv(file.path(RES,"main_preQC_cell_counts.csv"), show_col_types = FALSE) %>%
  mutate(sample = factor(sample, levels = SAMP))
p1 <- qc %>% select(sample, pre_QC_cells, post_QC_cells) %>%
  pivot_longer(-sample, names_to="stage", values_to="cells") %>%
  mutate(stage = recode(stage, pre_QC_cells="Pre-QC", post_QC_cells="Post-QC"),
         stage = factor(stage, levels=c("Pre-QC","Post-QC"))) %>%
  ggplot(aes(sample, cells, fill = sample, alpha = stage)) +
  geom_col(position = position_dodge(width=.75), width=.68, colour="grey30", linewidth=.2) +
  scale_fill_manual(values = sample_colors, guide="none") +
  scale_alpha_manual(values = c("Pre-QC"=.45, "Post-QC"=1), name=NULL,
                     guide = guide_legend(override.aes = list(fill = "#8C5A5A"))) +
  scale_y_continuous(labels = comma, expand = expansion(mult=c(0,.08))) +
  labs(title = "Cell numbers before and after QC",
       subtitle = "89,573 pre-QC  ->  76,366 post-QC (13,207 removed, 14.7%)",
       x = NULL, y = "Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_panel(p1, "P01_QC_cell_numbers", 4.3, 2.9)


at <- read_csv(file.path(RES,"main_raw_QC_threshold_attribution.csv"), show_col_types = FALSE)
names(at) <- c("sample","pre","post","feat_low","feat_high","umi_high","mt_high","hb_high")
THR <- c("percent.mt >= 10"="#C8956B", "nCount >= 80,000"="#B05A5A",
         "nFeature <= 200"="#6B8E9B", "nFeature >= 10,000"="#8BA6B2", "percent.HB >= 0.35"="#D4CFC5")
p2 <- at %>% mutate(sample = factor(sample, levels = SAMP)) %>%
  select(sample, `percent.mt >= 10`=mt_high, `nCount >= 80,000`=umi_high,
         `nFeature <= 200`=feat_low, `nFeature >= 10,000`=feat_high,
         `percent.HB >= 0.35`=hb_high) %>%
  pivot_longer(-sample, names_to="threshold", values_to="removed") %>%
  mutate(threshold = factor(threshold, levels = names(THR))) %>%
  ggplot(aes(sample, removed, fill = threshold)) +
  geom_col(width=.68, colour="grey30", linewidth=.2) +
  scale_fill_manual(values = THR, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult=c(0,.08))) +
  labs(title = "Cells removed by each QC threshold",
       subtitle = paste0("percent.mt accounts for 12,981 of 13,207 removals; ",
                         "the nFeature bounds remove 1 cell in total and percent.HB removes none"),
       x = NULL, y = "Cells removed") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.key.size = unit(.3,"cm"))
save_panel(p2, "P02_QC_threshold_attribution", 4.7, 2.9)


db <- read_csv(file.path(RES,"main_doublet_assessment.csv"), show_col_types = FALSE) %>%
  mutate(sample = factor(sample, levels = SAMP))
db2 <- db %>% select(sample, `score >= 0.20`=frac_ge_0.20, `score >= 0.33`=frac_ge_0.33,
                     `score >= 0.50`=frac_ge_0.50) %>%
  pivot_longer(-sample, names_to="cut", values_to="pct") %>%
  mutate(cut = factor(cut, levels=c("score >= 0.20","score >= 0.33","score >= 0.50")))
p3 <- ggplot() +
  geom_col(data = db2, aes(sample, pct, fill = sample, alpha = cut),
           position=position_dodge(width=.75), width=.66, colour="grey30", linewidth=.2) +
  geom_point(data = db, aes(sample, expected_rate_pct_10x,
                            shape = "10x Genomics expected\nmultiplet rate (interpolated)"),
             size=2.0, colour="grey15", stroke=.9, inherit.aes=FALSE) +
  scale_fill_manual(values = sample_colors, guide="none") +
  scale_alpha_manual(values=c(.45,.75,1), name="Artificial-doublet score",
                     guide = guide_legend(override.aes = list(fill = "#8C5A5A", shape = NA))) +
  scale_shape_manual(values = setNames(4, "10x Genomics expected\nmultiplet rate (interpolated)"),
                     name = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 2.4))) +
  scale_y_continuous(expand = expansion(mult=c(0,.08))) +
  labs(title = "Doublet assessment",
       subtitle = "Bars: share of cells above an artificial-doublet score cutoff.  Cross: 10x expected multiplet rate.",
       x = NULL, y = "Cells (%)") +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top",
        legend.box="vertical", legend.spacing.y=unit(0.05,"cm"))
save_panel(p3, "P03_doublet_assessment", 4.9, 3.2)





GRP_OF <- c("H-CD8Trp-1"="H-Trp","H-CD8Trp-2"="H-Trp",
            "E-CD8Trp-1"="E-Trp","E-CD8Trp-2"="E-Trp",
            "E-CD8DP-1"="E-DP","E-CD8DP-2"="E-DP")
ds <- read_csv(file.path(RES,"main_productive_TCR_equal_cell_summary.csv"), show_col_types=FALSE) %>%
  filter(definition == "CTstrict_like") %>%
  mutate(sample = factor(as.character(sample), levels=SAMP),
         donor  = factor(ifelse(grepl("-1$", as.character(sample)), "Donor 1", "Donor 2"),
                         levels=c("Donor 1","Donor 2")),
         group  = factor(unname(GRP_OF[as.character(sample)]), levels=GRP))
ds$group_num <- as.numeric(factor(ds$group, levels=GRP))
p4 <- ggplot(ds, aes(group_num, downsample_median)) +
  geom_line(aes(group = donor), colour="grey35", linewidth=.5, linetype="22") +
  geom_point(aes(fill = sample), shape=21, size=2.6, colour="grey20", stroke=.5) +
  geom_text(aes(label=sprintf("%.2f", downsample_median)), vjust=-1.4, size=1.9, colour="grey25") +
  facet_wrap(~ donor, nrow=1) +
  scale_fill_manual(values = sample_colors, name=NULL) +
  scale_x_continuous(breaks=1:3, labels=COND, limits=c(0.60, 3.40)) +
  scale_y_continuous(expand = expansion(mult = c(.10, .12))) +
  labs(title = "TCR Shannon diversity after equal-cell downsampling",
       subtitle = "n = 6,752 productive TCR-bearing cells per library; 2,000 draws without replacement",
       x = NULL, y = "Shannon diversity (median)") +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
save_panel(p4, "P04_TCR_equal_cell_downsampling", 5.4, 2.9)


dp <- read_csv(file.path(RES,"main_depth_matched_TCR_gene_level.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP), floor_umis=factor(floor_umis))
p5 <- ggplot(dp, aes(floor_umis, downsample_median, colour=sample, group=sample)) +
  geom_line(linewidth=.45) + geom_point(size=1.3) +
  scale_colour_manual(values = sample_colors, name=NULL) +
  labs(title = "Diversity across VDJ UMI floors",
       subtitle = "Floors 1 / 5 / 10 / 15 / 20 give n = 6,752 / 6,156 / 4,351 / 2,221 / 990",
       x = "Minimum VDJ UMIs per cell", y = "Shannon diversity (median)") +
  theme(legend.position="right")
save_panel(p5, "P05_depth_floor_sweep", 4.7, 2.9)


pb <- data.frame(
  comparison = factor(c("E-Trp vs H-Trp","E-DP vs E-Trp","E-DP vs H-Trp"),
                      levels=c("E-Trp vs H-Trp","E-DP vs E-Trp","E-DP vs H-Trp")),
  fdr=c(12748,15585,18193), lfc=c(9808,12241,14980))

PBCOL <- c("E-Trp vs H-Trp"="#B05A5A", "E-DP vs E-Trp"="#C8956B", "E-DP vs H-Trp"="#6B8E9B")
p6 <- pb %>% select(comparison, `FDR < 0.05`=fdr, `and |log2FC| > 1`=lfc) %>%
  pivot_longer(-comparison, names_to="criterion", values_to="n") %>%
  ggplot(aes(comparison, n, fill=comparison, alpha=criterion)) +
  geom_col(position=position_dodge(width=.72), width=.64, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=26083, linewidth=.3, linetype="dashed", colour="grey40") +
  annotate("text", x=.60, y=26083, label="26,083 genes tested", hjust=0, vjust=-.4, size=2) +
  scale_fill_manual(values=PBCOL, guide="none") +
  scale_alpha_manual(values=c(.45,1), name=NULL,
                     guide = guide_legend(override.aes = list(fill="#B05A5A"))) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0,.12))) +
  labs(title = "Donor-paired pseudobulk: significant genes",
       subtitle = "edgeR TMM + quasi-likelihood, ~ donor + group; two donors",
       x=NULL, y="Genes") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="top")
save_panel(p6, "P06_pseudobulk_significant_genes", 4.7, 2.9)


sg <- read_csv(file.path(RES,"main_recomputed_signature_donor_summary.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP),
         donor=factor(donor, levels=c("Donor1","Donor2")),
         group=factor(group, levels=GRP))
levels(sg$donor) <- DON
SIGORD <- c("memory4","exhaustion4","Exhaustion7","Stemness","TissueResident",
            "Glycolysis","OXPHOS","Costimulation","Migration",
            "Pathogenic_cytotoxic","CD8_Treg_associated","CD8_Trp_50")
sg$signature <- factor(sg$signature, levels = SIGORD)
p7 <- ggplot(sg, aes(group, mean_control_adjusted_score, fill=sample)) +
  geom_col(width=.7, colour="grey30", linewidth=.18) +
  geom_hline(yintercept=0, linewidth=.25, colour="grey30") +
  facet_wrap(~ signature, ncol=4, scales="free_y") +
  scale_fill_manual(values=sample_colors, name=NULL) +
  scale_x_discrete(labels=COND) +
  labs(title="Signature scores at donor level",
       subtitle="Control-adjusted mean per donor-condition; no cell-level p-values are reported",
       x=NULL, y="Control-adjusted mean score") +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=5.5),
        strip.text=element_text(size=6), legend.position="bottom")
save_panel(p7, "P07_signature_donor_level", 6.6, 4.4)


eg <- read_csv(file.path(RES,"main_requested_eight_genes_donor_expression.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP))
ord <- c("TCF7","CCR7","IL7R","LEF1","TOX","HAVCR2","ENTPD1","LAG3")
eg$gene <- factor(eg$gene, levels=rev(ord))
p8 <- ggplot(eg, aes(sample, gene, fill=mean_log1p_normalized_RNA)) +
  geom_tile(colour="white", linewidth=.4) +
  geom_text(aes(label=sprintf("%.2f", mean_log1p_normalized_RNA)), size=1.6, colour="grey15") +
  scale_fill_gradientn(colours=c("#D4CFC5", "#B5C1B4", "#8C5A5A"), name="log1p") +
  labs(title="Memory and exhaustion genes",
       subtitle="Mean log1p normalised RNA per donor-condition", x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())
save_panel(p8, "P08_memory_exhaustion_genes", 4.6, 2.9)


sc <- read_csv(file.path(RES,"main_requested_memory_exhaustion_donor_scores.csv"), show_col_types=FALSE) %>%
  mutate(sample=factor(sample, levels=SAMP), donor=factor(donor, levels=c("Donor1","Donor2")),
         group=factor(group, levels=GRP))
levels(sc$donor) <- DON
p9 <- sc %>% select(sample, group, donor,
                    `Memory (4-gene)`=memory4_control_adjusted_mean,
                    `Exhaustion (4-gene)`=exhaustion4_control_adjusted_mean) %>%
  pivot_longer(c(`Memory (4-gene)`, `Exhaustion (4-gene)`), names_to="score", values_to="value") %>%
  ggplot(aes(group, value, fill=sample)) +
  geom_col(width=.66, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey30") +
  facet_wrap(~ score + donor, nrow = 1, scales="free_y") +
  scale_fill_manual(values=sample_colors, name=NULL) +
  scale_x_discrete(labels=COND) +
  labs(title="Memory and exhaustion scores at donor level",
       subtitle="Control-adjusted, AddModuleScore (nbin 24, ctrl 100, seed 20261004)",
       x=NULL, y="Control-adjusted mean score") +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
save_panel(p9, "P09_memory_exhaustion_scores", 5.2, 3.1)


tr <- read_csv(file.path(RES,"main_productive_TCR_donor_trends.csv"), show_col_types=FALSE) %>%
  mutate(donor = factor(donor, levels=c("Donor1","Donor2")),
         comparison = factor(comparison, levels=c("Pre-sel vs SS","Post-sel vs Pre-sel")),
         definition = recode(definition, CTgene_like="CTgene", CTstrict_like="CTstrict",
                             .default = definition))
levels(tr$donor) <- DON
p10 <- ggplot(tr, aes(comparison, downsample_delta, fill=donor)) +
  geom_col(position=position_dodge(width=.72), width=.64, colour="grey30", linewidth=.2) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey30") +
  facet_wrap(~ definition, nrow=1, scales="free_y") +
  scale_fill_manual(values=c("Donor 1"="#8C5A5A","Donor 2"="#A87575"), name=NULL) +
  labs(title="Direction of change by donor",
       subtitle="All 8 of 8 adjacent comparisons keep the sign of the raw result",
       x=NULL, y="Change in downsampled Shannon diversity") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="top")
save_panel(p10, "P10_per_donor_direction", 5.2, 2.9)

cat("\nEnglish panels written to", OUT, ":", length(list.files(OUT, pattern="\\.pdf$")), "PDF\n")
