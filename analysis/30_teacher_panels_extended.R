








suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
setwd(Sys.getenv("DATA_ROOT", "."))
OUT <- "figures/teacher_panels"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

PANELS <- list(
 "Exhaustion" = c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT",
                  "TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1"),
 "Activation / effector function" = c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1",
                  "KLRD1","KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2",
                  "ID2","ID3","PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1",
                  "NR4A2","NR4A3","PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB",
                  "JUND","STAT1","STAT2","STAT5A","STAT6","STAT4","EOMES"),
 "Cytotoxicity" = c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1","IFNG","TNF","SERPINB1",
                  "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3",
                  "CST7","CSTB","LAMP1","LAMP3","CAPN2"))
GRP  <- c("H-Trp","E-Trp","E-DP"); COND <- c("SS","Pre-sel","Post-sel")

PANEL_LEVELS <- c("Exhaustion","Activation / effector function","Cytotoxicity")
DONCOL <- c("Donor 1"="#C8956B", "Donor 2"="#6B8E9B")
UPDOWN <- c("up in both"="#B05A5A", "down in both"="#6B8E9B", "discordant"="#B9B4AC")
YEL <- c("#FFFDE7","#FFF9C4","#FFE082","#FFC107","#FF8F00")

o <- readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
o$group <- factor(o$group, levels=GRP)
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
cat("  :", ncol(o), "\n")


for (nm in names(PANELS)) {
  g <- intersect(PANELS[[nm]], rownames(o))
  cat("", nm, ": ", length(g), "/", length(PANELS[[nm]]), "\n")
  o <- AddModuleScore(o, features=list(g), name=paste0("SC_", gsub("[^A-Za-z]", "", nm)),
                      nbin=24, ctrl=100, seed=20261004)
}
score_cols <- grep("^SC_", colnames(o@meta.data), value=TRUE)

stopifnot(length(score_cols) == length(PANEL_LEVELS))
lab <- setNames(PANEL_LEVELS, score_cols)

sc <- FetchData(o, vars=c(score_cols, "group", "donor"), layer="data") %>%
  pivot_longer(all_of(score_cols), names_to="panel_col", values_to="score") %>%
  mutate(panel=factor(unname(lab[panel_col]), levels=PANEL_LEVELS)) %>%
  group_by(donor, group, panel) %>%
  summarise(mean_score=mean(score), sd_score=sd(score), n=n(), .groups="drop")
write.csv(sc, "results/teacher_panel_module_scores.csv", row.names=FALSE)


pA <- ggplot(sc, aes(group, mean_score, colour=donor, group=donor)) +
  geom_line(linewidth=.5) + geom_point(size=2.2) +
  facet_wrap(~ panel, scales="free_y", nrow=1) +
  scale_colour_manual(values=DONCOL, name=NULL) +
  labs(title="Panel module scores across the three conditions",
       subtitle="AddModuleScore per panel (nbin 24, ctrl 100, seed 20261004); donor-level means",
       x=NULL, y="mean module score") +
  theme_bw(base_size=9) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
ggsave(file.path(OUT, "G7_panel_module_scores.pdf"), pA, width=8.4, height=3.2, device=cairo_pdf)
cat("  : G7_panel_module_scores.pdf\n")


d <- read.csv("results/teacher_panel_donor_delta.csv", fileEncoding="UTF-8-BOM")
d$gene <- as.character(d$gene)
d$panel <- factor(d$panel, levels=PANEL_LEVELS)
dd <- d %>% select(gene, panel, d1_post=p1, d2_post=p2, mean_post) %>%
  pivot_longer(c(d1_post, d2_post), names_to="who", values_to="delta") %>%
  mutate(donor=ifelse(who=="d1_post","Donor 1","Donor 2"))
pB <- ggplot(dd, aes(delta, reorder(gene, mean_post), colour=donor)) +
  geom_vline(xintercept=0, colour="grey40", linetype="22", linewidth=.35) +
  geom_point(size=1.7) +
  facet_wrap(~ panel, scales="free_y", nrow=1) +
  scale_colour_manual(values=DONCOL, name=NULL) +
  labs(title="Panel markers: Post-sel against Pre-sel, donor level",
       subtitle="Δ mean log1p (E-DP − E-Trp); one point per donor",
       x="Δ mean log1p (Post-sel − Pre-sel)", y=NULL) +
  theme_bw(base_size=8) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        axis.text.y=element_text(size=5.6),
        plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
ggsave(file.path(OUT, "G8_post_vs_pre_delta.pdf"), pB, width=9.4, height=max(4.0, .13*nrow(d)*3/3+2.4), device=cairo_pdf)
cat("  : G8_post_vs_pre_delta.pdf\n")


emb <- Embeddings(o, "umap")
key <- d %>% group_by(panel) %>% slice_max(abs(mean_delta), n=6) %>% ungroup()
expr <- FetchData(o, vars=key$gene, layer="data")
for (p in names(PANELS)) {
  gs <- key %>% filter(panel==p) %>% pull(gene)
  pl <- lapply(gs, function(g) {
    v <- expr[[g]]
    lo <- quantile(v, .05); hi <- quantile(v, .95)
    df <- data.frame(x=emb[,1], y=emb[,2], v=pmin(pmax(v, lo), hi))
    ggplot(df, aes(x, y, colour=v)) + geom_point(size=.06) +
      scale_colour_gradientn(colours=YEL, name=NULL) +
      labs(title=g, x=NULL, y=NULL) +
      theme_void(base_size=7) +
      theme(plot.title=element_text(size=7.5, hjust=.5), legend.key.width=unit(.28,"cm"),
            legend.key.height=unit(.28,"cm"), legend.text=element_text(size=5))
  })
  fname <- sprintf("G9%d_FeaturePlot_%s.pdf", which(names(PANELS)==p),
                   c("exhaustion","activation","cytotoxicity")[which(names(PANELS)==p)])
  g <- wrap_plots(pl, ncol=6) +
    plot_annotation(title=sprintf("%s: key markers on the UMAP", p),
                    subtitle="q5-q95 colour scale, point size 0.06; all 76,366 QC-passed cells",
                    theme=theme(plot.title=element_text(face="bold", size=10),
                                plot.subtitle=element_text(size=8)))
  ggsave(file.path(OUT, fname), g, width=9.2, height=1.9, device=cairo_pdf)

  ggsave(file.path(OUT, sub("\\.pdf$", ".png", fname)), g, width=9.2, height=1.9, dpi=600)
  cat("  :", fname, "+ PNG\n")
}
cat("\n  。:", OUT, "\n")
