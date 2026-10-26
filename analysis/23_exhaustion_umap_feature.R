








suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
  library(ggrastr)
})
setwd("DATA_ROOT")
OUT <- "figures/exhaustion"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


morandi_group <- c("E-DP"="#8C5A5A", "E-Trp"="#6B7B6B", "H-Trp"="#7A8B99")
group_colors  <- c("E-DP"="#C8956B", "E-Trp"="#B05A5A", "H-Trp"="#6B8E9B")
sample_colors <- c("E-CD8DP-1"="#C8956B", "E-CD8DP-2"="#D9AE8A",
                   "E-CD8Trp-1"="#B05A5A", "E-CD8Trp-2"="#C27A7A",
                   "H-CD8Trp-1"="#6B8E9B", "H-CD8Trp-2"="#8BA6B2")
RAMP <- c("#D4CFC5", "#B5C1B4", "#8C5A5A")
COND <- c("H-Trp"="SS", "E-Trp"="Pre-sel", "E-DP"="Post-sel")
SAMP <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")

theme_pub <- theme_minimal(base_size = 7) +
  theme(plot.title = element_text(size = 8, face="bold", hjust=.5, margin=margin(b=3)),
        plot.subtitle = element_text(size = 7, hjust=.5),
        axis.title = element_text(size=7), axis.text = element_text(size=6, color="black"),
        axis.line = element_line(linewidth=.3, color="black"),
        axis.ticks = element_line(linewidth=.3, color="black"),
        legend.title = element_text(size=7), legend.text = element_text(size=6),
        legend.key.size = unit(.4,"cm"),
        panel.grid.major = element_line(linewidth=.25, color="grey90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth=.3, fill=NA, color="black"))
theme_set(theme_pub)
save_panel <- function(p, name, w=4.2, h=3.1) {
  cairo_pdf(file.path(OUT, paste0(name,".pdf")), width=w, height=h); print(p); dev.off()
  ggsave(file.path(OUT, paste0(name,".png")), p, width=w, height=h, dpi=600)
  cat("  :", name, "\n")
}


cat(" seurat_filtered_newQC.rds …\n")
o <- readRDS("DATA_ROOT")
o$sample <- factor(o$sample, levels = SAMP)
o$group  <- factor(o$group,  levels = c("H-Trp","E-Trp","E-DP"))
o$donor  <- factor(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"),
                   levels = c("Donor 1","Donor 2"))
o$condition <- factor(unname(COND[as.character(o$group)]), levels = c("SS","Pre-sel","Post-sel"))


MEM <- c("TCF7","CCR7","IL7R","LEF1"); EXH <- c("TOX","HAVCR2","ENTPD1","LAG3")
set.seed(20261004)
o <- AddModuleScore(o, list(intersect(MEM, rownames(o))), name="memory4", nbin=24, ctrl=100, search=FALSE)
set.seed(20261004)
o <- AddModuleScore(o, list(intersect(EXH, rownames(o))), name="exhaustion4", nbin=24, ctrl=100, search=FALSE)
o$memory4 <- o[["memory41"]]; o[["memory41"]] <- NULL
o$exhaustion4 <- o[["exhaustion41"]]; o[["exhaustion41"]] <- NULL
cat("  : memory4 / exhaustion4（AddModuleScore  *1，）\n")


p <- DimPlot(o, group.by="sample", reduction="umap", cols=sample_colors, pt.size=.25) +
  labs(title="UMAP by sample") + theme(legend.position="right", legend.key.size=unit(.28,"cm"))
p <- ggrastr::rasterise(p, layers="GeomPoint", dpi=300)
save_panel(p, "U01_UMAP_by_sample", 4.6, 3.0)

p <- DimPlot(o, group.by="condition", reduction="umap", cols=setNames(group_colors, c("SS","Pre-sel","Post-sel")),
             pt.size=.25) + labs(title="UMAP by condition") +
  theme(legend.position="right")
p <- ggrastr::rasterise(p, layers="GeomPoint", dpi=300)
save_panel(p, "U02_UMAP_by_condition", 4.4, 3.0)

p <- DimPlot(o, group.by="donor", reduction="umap",
             cols=c("Donor 1"="#8C5A5A","Donor 2"="#A87575"), pt.size=.25) +
  labs(title="UMAP by donor") + theme(legend.position="right")
p <- ggrastr::rasterise(p, layers="GeomPoint", dpi=300)
save_panel(p, "U03_UMAP_by_donor", 4.2, 3.0)

p <- DimPlot(o, group.by="seurat_clusters", reduction="umap", label=TRUE, repel=TRUE, pt.size=.25) +
  labs(title="UMAP by cluster") + theme(legend.position="none")
p <- ggrastr::rasterise(p, layers="GeomPoint", dpi=300)
save_panel(p, "U04_UMAP_by_cluster", 4.2, 3.0)


fp <- function(genes, name, title, ncol=2, w=4.6, h=4.4) {
  g <- intersect(genes, rownames(o))
  pl <- FeaturePlot(o, features=g, reduction="umap", min.cutoff="q5", max.cutoff="q95",
                    pt.size=.5, combine=FALSE)
  pl <- lapply(pl, function(x) ggrastr::rasterise(
    x + scale_color_gradientn(colors=RAMP, name="expr") + theme_pub +
      theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()),
    layers="GeomPoint", dpi=300))
  p <- wrap_plots(pl, ncol=ncol)
  save_panel(p, name, w, h)
}
fp(MEM, "F01_FeaturePlot_memory",    "Memory / progenitor genes")
fp(EXH, "F02_FeaturePlot_exhaustion","Exhaustion-associated genes")


emb <- as.data.frame(Embeddings(o, "umap")); colnames(emb) <- c("UMAP1","UMAP2")
emb$memory4 <- o$memory4; emb$exhaustion4 <- o$exhaustion4
sc <- do.call(rbind, lapply(c("memory4","exhaustion4"), function(k)
  data.frame(emb[,c("UMAP1","UMAP2")], score=emb[[k]], signature=k)))
sc$signature <- factor(sc$signature, levels=c("memory4","exhaustion4"))
p <- ggplot(sc, aes(UMAP1, UMAP2, colour=score)) +
  geom_point(size=.22) +
  facet_wrap(~ signature, nrow=1) +
  scale_colour_gradientn(colours=RAMP, name="score",
                         limits=quantile(sc$score, c(.05,.95)), oob=scales::squish) +
  labs(title="Module scores on UMAP") + theme_pub +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
        strip.text=element_text(size=7, face="bold"))
p <- ggrastr::rasterise(p, layers="GeomPoint", dpi=300)
save_panel(p, "F03_FeaturePlot_scores", 5.0, 2.9)


vp <- function(features, name, title, cols, ncol=2, w=5.0, h=4.4) {
  p <- VlnPlot(o, features=features, group.by="group", pt.size=0, cols=cols,
               ncol=ncol, combine=TRUE) & theme_pub &
    theme(axis.title.x=element_blank(), legend.position="none")
  p <- p + plot_annotation(title=title, theme=theme(plot.title=element_text(size=9, face="bold", hjust=.5)))
  save_panel(p, name, w, h)
}
vp(MEM, "V01_violin_memory",     "Memory / progenitor genes",     group_colors, 2, 5.0, 4.4)
vp(EXH, "V02_violin_exhaustion", "Exhaustion-associated genes",   group_colors, 2, 5.0, 4.4)
vp(c("memory4","exhaustion4"), "V03_violin_scores", "Module scores by condition", group_colors, 2, 4.4, 3.4)
cat("\n，:", OUT, "\n")
