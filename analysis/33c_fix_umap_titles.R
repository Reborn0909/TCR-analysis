

suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork)})
setwd("DATA_ROOT")
OUT <- "figures/teacher_panels_umap_dot"
PANELS <- list(
 "Exhaustion" = c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT","TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1"),
 "Activation / effector function" = c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1","KLRD1","KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2","ID2","ID3","PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1","NR4A2","NR4A3","PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB","JUND","STAT1","STAT2","STAT5A","STAT6","STAT4","EOMES"),
 "Cytotoxicity" = c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1","IFNG","TNF","SERPINB1","SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7","CSTB","LAMP1","LAMP3","CAPN2"))
GRP <- c("H-Trp","E-Trp","E-DP")
YELLOW <- c("#FFFDE7","#FFF9C4","#FFE082","#FFC107","#FF8F00")
tag <- c("Exhaustion"="exhaustion","Activation / effector function"="activation_effector","Cytotoxicity"="cytotoxicity")

o <- readRDS("DATA_ROOT")
o$group <- factor(o$group, levels=GRP)
emb  <- Embeddings(o, "umap")
expr <- FetchData(o, vars=unique(unlist(PANELS)), layer="data")


q95 <- sapply(colnames(expr), function(g) as.numeric(quantile(expr[[g]], .95)))
cat(sprintf("  q95  %.3f ~ %.2f\n", min(q95), max(q95)))

for (pn in names(PANELS)) {
  gs <- intersect(PANELS[[pn]], rownames(o))
  pl <- lapply(gs, function(g) {
    v <- expr[[g]]; lo <- quantile(v,.05); hi <- q95[[g]]
    df <- data.frame(x=emb[,1], y=emb[,2], v=pmin(pmax(v,lo),hi))
    ggplot(df, aes(x,y,colour=v)) + geom_point(size=.05) +
      scale_colour_gradientn(colours=YELLOW, name=NULL, guide="none") +
      labs(title=sprintf("%s  (q95 %.2f)", g, hi)) +
      theme_void(base_size=7) +
      theme(plot.title=element_text(size=7, hjust=.5))
  })
  ncol <- min(7, length(pl)); nrow <- ceiling(length(pl)/ncol)
  gg <- wrap_plots(pl, ncol=ncol) +
    plot_annotation(title=sprintf("%s: all %d markers on the UMAP", pn, length(gs)),
                    subtitle="light-to-deep yellow, clipped at the 5th and 95th percentiles; the q95 value in each title is that gene's 95th percentile, so a pale panel means low expression, not uniform expression; 76,366 QC-passed cells",
                    theme=theme(plot.title=element_text(face="bold", size=10),
                                plot.subtitle=element_text(size=8)))
  f <- sprintf("U%d_UMAP_all_%s.pdf", which(names(PANELS)==pn), tag[[pn]])
  ggsave(file.path(OUT, f), gg, width=1.32*ncol, height=1.36*nrow, device=cairo_pdf)
  ggsave(file.path(OUT, sub("\\.pdf$", ".png", f)), gg, width=1.32*ncol, height=1.36*nrow, dpi=600)
  cat("  :", f, "+ PNG\n")
}
cat("  \n")
