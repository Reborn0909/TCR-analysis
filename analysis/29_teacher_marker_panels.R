







suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
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
ALLG <- unique(unlist(PANELS))
cat("  :", length(ALLG), "\n")


DONCOL <- c("Donor 1"="#C8956B", "Donor 2"="#6B8E9B")
UPDOWN <- c("up in both"="#B05A5A", "down in both"="#6B8E9B", "discordant"="#B9B4AC")
GRP    <- c("H-Trp","E-Trp","E-DP")
COND   <- c("SS","Pre-sel","Post-sel")
HEAT   <- c("#D4CFC5","#B5C1B4","#8C5A5A")

o <- readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
o$group <- factor(o$group, levels=GRP)
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
cat("  :", ncol(o), " ·", nrow(o), "\n")

genes <- intersect(ALLG, rownames(o))
miss  <- setdiff(ALLG, genes)
cat("  :", length(genes), "/", length(ALLG), "\n")
if (length(miss)) cat("  **:", paste(miss, collapse=", "), "\n")


emb <- FetchData(o, vars=c(genes, "group", "donor"), layer="data")
long <- emb %>% pivot_longer(all_of(genes), names_to="gene", values_to="expr") %>%
  group_by(donor, group, gene) %>%
  summarise(mean_log1p = mean(expr), pct_expressing = 100*mean(expr>0),
            n_cells = n(), .groups="drop")
write.csv(long, "results/teacher_panel_donor_expression.csv", row.names=FALSE)

w <- long %>% select(donor, group, gene, mean_log1p) %>%
  pivot_wider(names_from=group, values_from=mean_log1p)
d <- w %>% mutate(delta_Pre_vs_SS  = `E-Trp` - `H-Trp`,
                  delta_Post_vs_Pre = `E-DP` - `E-Trp`) %>%
  select(donor, gene, delta_Pre_vs_SS, delta_Post_vs_Pre)

agg <- d %>% group_by(gene) %>%
  summarise(d1 = delta_Pre_vs_SS[donor=="Donor 1"], d2 = delta_Pre_vs_SS[donor=="Donor 2"],
            mean_delta = mean(delta_Pre_vs_SS),
            p1 = delta_Post_vs_Pre[donor=="Donor 1"], p2 = delta_Post_vs_Pre[donor=="Donor 2"],
            mean_post = mean(delta_Post_vs_Pre),
            up_both   = all(delta_Pre_vs_SS > 0),
            down_both = all(delta_Pre_vs_SS < 0),
            .groups="drop") %>%
  mutate(panel = sapply(gene, function(g) names(PANELS)[sapply(PANELS, function(x) g %in% x)][1]),
         call  = ifelse(up_both, "up in both",
                 ifelse(down_both, "down in both", "discordant"))) %>%
  arrange(panel, desc(mean_delta))
write.csv(agg, "results/teacher_panel_donor_delta.csv", row.names=FALSE)


sm <- agg %>% group_by(panel) %>%
  summarise(n = n(), up = sum(up_both), down = sum(down_both),
            disc = sum(!up_both & !down_both), .groups="drop")
cat("\n===  （Pre-sel vs SS）===\n")
for (i in seq_len(nrow(sm))) cat(sprintf("  %-32s n=%2d   %2d   %2d   %d\n",
  sm$panel[i], sm$n[i], sm$up[i], sm$down[i], sm$disc[i]))

cat("\n=== （）===\n")
for (p in names(PANELS)) {
  u <- agg %>% filter(panel==p, up_both) %>% arrange(desc(mean_delta))
  if (nrow(u)) for (i in seq_len(nrow(u)))
    cat(sprintf("  %-32s %-10s D1 %+0.3f  D2 %+0.3f   %+0.3f\n",
        p, u$gene[i], u$d1[i], u$d2[i], u$mean_delta[i]))
}
cat("\n===  ===\n")
for (p in names(PANELS)) {
  u <- agg %>% filter(panel==p, down_both) %>% arrange(mean_delta)
  if (nrow(u)) for (i in seq_len(nrow(u)))
    cat(sprintf("  %-32s %-10s D1 %+0.3f  D2 %+0.3f   %+0.3f\n",
        p, u$gene[i], u$d1[i], u$d2[i], u$mean_delta[i]))
}
sac <- agg %>% filter(!up_both & !down_both)
if (nrow(sac)) { cat("\n===  ===\n")
  for (i in seq_len(nrow(sac))) cat(sprintf("  %-32s %-10s D1 %+0.3f  D2 %+0.3f\n",
    sac$panel[i], sac$gene[i], sac$d1[i], sac$d2[i])) }


mk_delta <- function(p, fname, w=6.6) {
  a <- agg %>% filter(panel==p) %>%
    mutate(gene=factor(gene, levels=gene[order(mean_delta)]))
  dd <- d %>% filter(gene %in% a$gene)
  p1 <- ggplot() +
    geom_segment(data=a, aes(y=gene, yend=gene, x=d1, xend=d2),
                 colour="grey60", linewidth=.45) +
    geom_point(data=dd, aes(delta_Pre_vs_SS, gene, fill=donor),
               shape=21, size=2.5, colour="grey25", stroke=.4) +
    geom_vline(xintercept=0, colour="grey35", linewidth=.35, linetype="22") +
    scale_fill_manual(values=DONCOL, name=NULL) +
    labs(title=sprintf("%s markers: Pre-sel vs SS, donor level", p),
         subtitle="Δ mean log1p normalised expression (E-Trp − H-Trp); a line joins the two donors",
         x="Δ mean log1p (Pre-sel − SS)", y=NULL) +
    theme_bw(base_size=9) +
    theme(panel.grid.minor=element_blank(), legend.position="top",
          plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
  ggsave(file.path(OUT, fname), p1, width=w, height=max(2.6, .22*nrow(a)+1.5), device=cairo_pdf)
  cat("  :", fname, "\n")
}
mk_delta("Exhaustion", "G1_exhaustion_delta.pdf", 6.2)
mk_delta("Activation / effector function", "G2_activation_effector_delta.pdf", 6.8)
mk_delta("Cytotoxicity", "G3_cytotoxicity_delta.pdf", 6.4)


hm <- long %>% mutate(col = paste(donor, group, sep=" · ")) %>%
  select(gene, col, mean_log1p)
ord <- agg %>% arrange(panel, desc(mean_delta)) %>% pull(gene)
hm$gene <- factor(hm$gene, levels=rev(ord))


hm$col  <- factor(hm$col, levels=as.vector(t(outer(c("Donor 1","Donor 2"), GRP, paste, sep=" · "))))
hm <- hm %>% group_by(gene) %>% mutate(z = as.numeric(scale(mean_log1p))) %>% ungroup() %>%
  mutate(z = ifelse(is.finite(z), z, 0))
p4 <- ggplot(hm, aes(col, gene, fill=z)) + geom_tile(colour="white", linewidth=.25) +
  scale_fill_gradientn(colours=HEAT, limits=c(-2.4,2.4), oob=scales::squish, name="row z") +
  labs(title="Teacher panels: donor-level mean expression",
       subtitle="76 markers · row z-score of mean log1p · columns are donor × condition",
       x=NULL, y=NULL) +
  theme_minimal(base_size=6.4) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6.4), axis.text.y=element_text(size=5.2),
        panel.grid=element_blank(), plot.title=element_text(face="bold", size=10),
        plot.subtitle=element_text(size=8))
ggsave(file.path(OUT, "G4_all_panels_heatmap.pdf"), p4, width=4.6, height=11.5, device=cairo_pdf)
cat("  : G4_all_panels_heatmap.pdf\n")


p5 <- ggplot(sm %>% pivot_longer(c(up,down,disc), names_to="k", values_to="v") %>%
               mutate(k=factor(k, levels=c("up","down","disc"),
                               labels=c("up in both","down in both","discordant"))),
             aes(panel, v, fill=k)) +
  geom_col(position="dodge", width=.7, colour="grey25", linewidth=.3) +
  geom_text(aes(label=v), position=position_dodge(width=.7), vjust=-.5, size=2.6) +
  scale_fill_manual(values=UPDOWN, name=NULL) +
  labs(title="Marker panels: how the two donors agree (Pre-sel vs SS)",
       subtitle="A marker counts as up or down only when both donors move the same way",
       x=NULL, y="Markers") +
  theme_bw(base_size=9) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        axis.text.x=element_text(angle=18,hjust=1),
        plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
ggsave(file.path(OUT, "G5_panel_summary.pdf"), p5, width=6.4, height=3.4, device=cairo_pdf)
cat("  : G5_panel_summary.pdf\n")


top <- agg %>% group_by(panel) %>% slice_max(abs(mean_delta), n=12) %>% ungroup()
tr <- long %>% filter(gene %in% top$gene) %>%
  mutate(gene=factor(gene, levels=top$gene), cond=as.numeric(group)) %>%
  left_join(top %>% select(gene, panel), by="gene")
p6 <- ggplot(tr, aes(cond, mean_log1p, colour=donor, group=donor)) +
  geom_line(linewidth=.45) + geom_point(size=1.5) +
  facet_wrap(~ gene, scales="free_y", ncol=6) +
  scale_colour_manual(values=DONCOL, name=NULL) +
  scale_x_continuous(breaks=1:3, labels=COND) +
  labs(title="Teacher panels: SS → Pre-sel → Post-sel, donor level",
       subtitle="The twelve largest movers in each panel", x=NULL, y="mean log1p") +
  theme_bw(base_size=8) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        axis.text.x=element_text(angle=30,hjust=1,size=6.4),
        plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
ggsave(file.path(OUT, "G6_trajectories_top_movers.pdf"), p6, width=9.2, height=6.4, device=cairo_pdf)
cat("  : G6_trajectories_top_movers.pdf\n")
cat("\n  。:", OUT, "\n")
