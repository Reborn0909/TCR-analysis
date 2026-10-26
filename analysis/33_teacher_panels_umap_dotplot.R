








suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
setwd(Sys.getenv("DATA_ROOT", "."))
OUT <- "figures/teacher_panels_umap_dot"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

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
PANEL_LEVELS <- names(PANELS)


GROUP_COLORS <- c("H-Trp"="#6B8E9B", "E-Trp"="#B05A5A", "E-DP"="#C8956B")
DON_COLORS   <- c("Donor 1"="#C8956B", "Donor 2"="#6B8E9B")
MORANDI_RAMP <- c("#D4CFC5","#B5C1B4","#8C5A5A")
YELLOW_RAMP  <- c("#FFFDE7","#FFF9C4","#FFE082","#FFC107","#FF8F00")

o <- readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
o$group <- factor(o$group, levels=GRP)
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
o$dc    <- factor(paste(o$donor, as.character(o$group), sep=" / "),
                  levels=as.vector(t(outer(c("Donor 1","Donor 2"), GRP, paste, sep=" / "))))
cat("  :", ncol(o), " · :", nlevels(o$dc), "\n")

tag <- c("Exhaustion"="exhaustion","Activation / effector function"="activation_effector",
         "Cytotoxicity"="cytotoxicity")


for (pn in PANEL_LEVELS) {
  g <- intersect(PANELS[[pn]], rownames(o))
  p <- DotPlot(o, features=g, group.by="dc", cols=MORANDI_RAMP, scale=TRUE,
               dot.scale=4.2) +
    RotatedAxis() +
    labs(title=sprintf("%s: expression across donor and condition", pn),
         subtitle="dot size = fraction of cells expressing; colour = scaled mean expression",
         x=NULL, y=NULL) +
    theme_bw(base_size=9) +
    theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=45,hjust=1,size=7),
          plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
  f <- sprintf("D%d_DotPlot_%s.pdf", which(PANEL_LEVELS==pn), tag[[pn]])
  ggsave(file.path(OUT, f), p, width=max(4.4, .24*length(g)+1.6), height=3.1, device=cairo_pdf)
  cat("  :", f, "\n")
}


emb  <- Embeddings(o, "umap")
expr <- FetchData(o, vars=unique(unlist(PANELS)), layer="data")
for (pn in PANEL_LEVELS) {
  gs <- intersect(PANELS[[pn]], rownames(o))
  pl <- lapply(gs, function(g) {
    v <- expr[[g]]; lo <- quantile(v,.05); hi <- quantile(v,.95)
    df <- data.frame(x=emb[,1], y=emb[,2], v=pmin(pmax(v,lo),hi))
    ggplot(df, aes(x,y,colour=v)) + geom_point(size=.05) +
      scale_colour_gradientn(colours=YELLOW_RAMP, name=NULL, guide="none") +
      labs(title=g) + theme_void(base_size=7) +
      theme(plot.title=element_text(size=7, hjust=.5))
  })
  ncol <- min(7, length(pl)); nrow <- ceiling(length(pl)/ncol)
  g <- wrap_plots(pl, ncol=ncol) +
    plot_annotation(title=sprintf("%s: all %d markers on the UMAP", pn, length(gs)),
                    subtitle="light-to-deep yellow, clipped at the 5th and 95th percentiles; 76,366 QC-passed cells",
                    theme=theme(plot.title=element_text(face="bold", size=10),
                                plot.subtitle=element_text(size=8)))
  f <- sprintf("U%d_UMAP_all_%s.pdf", which(PANEL_LEVELS==pn), tag[[pn]])
  ggsave(file.path(OUT, f), g, width=1.32*ncol, height=1.32*nrow, device=cairo_pdf)
  ggsave(file.path(OUT, sub("\\.pdf$", ".png", f)), g, width=1.32*ncol, height=1.32*nrow, dpi=600)
  cat("  :", f, "+ PNG（", length(gs), " ）\n")
}


for (pn in PANEL_LEVELS) {
  g <- intersect(PANELS[[pn]], rownames(o))
  col <- paste0("SC_", gsub("[^A-Za-z]", "", pn))
  if (!col %in% colnames(o@meta.data)) {
    o <- AddModuleScore(o, features=list(g), name=col, nbin=24, ctrl=100, seed=20261004)
    col <- paste0(col, "1")
  }
  df <- FetchData(o, vars=c(col, "group", "donor"), layer="data")
  names(df)[1] <- "score"
  p <- ggplot(df, aes(group, score, fill=group)) +
    geom_violin(scale="width", linewidth=.25, colour="grey30") +
    geom_boxplot(width=.14, outlier.shape=NA, linewidth=.25, colour="grey25") +
    facet_wrap(~ donor, nrow=1) +
    scale_fill_manual(values=GROUP_COLORS, name=NULL) +
    scale_x_discrete(labels=COND) +
    labs(title=sprintf("%s module score by condition and donor", pn),
         subtitle="AddModuleScore (nbin 24, ctrl 100, seed 20261004); cells shown for display only, no cell-level test",
         x=NULL, y="module score") +
    theme_bw(base_size=9) +
    theme(panel.grid.minor=element_blank(), legend.position="top",
          plot.title=element_text(face="bold", size=10), plot.subtitle=element_text(size=8))
  f <- sprintf("V%d_Vln_score_%s.pdf", which(PANEL_LEVELS==pn), tag[[pn]])
  ggsave(file.path(OUT, f), p, width=6.6, height=3.4, device=cairo_pdf)
  cat("  :", f, "\n")
}
cat("\n  。:", OUT, "\n")
