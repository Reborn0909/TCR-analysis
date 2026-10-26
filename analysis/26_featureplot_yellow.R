





suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
setwd("DATA_ROOT")
OUT <- "figures/featureplot_yellow"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


YEL <- c("#FFFDE7", "#FFF9C4", "#FFE082", "#FFC107", "#FF8F00")

theme_pub <- theme_minimal(base_size=7) +
  theme(plot.title=element_text(size=8.5, face="bold", hjust=.5, margin=margin(b=3)),
        plot.subtitle=element_text(size=7, hjust=.5),
        axis.title=element_text(size=7), axis.text=element_text(size=6, color="black"),
        axis.line=element_line(linewidth=.3, color="black"),
        axis.ticks=element_line(linewidth=.3, color="black"),
        legend.title=element_text(size=7), legend.text=element_text(size=6),
        legend.key.size=unit(.4,"cm"),
        panel.grid.major=element_line(linewidth=.25, color="grey90"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(linewidth=.3, fill=NA, color="black"),
        strip.text=element_text(size=7, face="bold"),
        strip.background=element_rect(linewidth=.3, fill="grey95", color="black"))
theme_set(theme_pub)

save_panel <- function(p, name, w, h) {
  cairo_pdf(file.path(OUT, paste0(name,".pdf")), width=w, height=h); print(p); dev.off()
  ggsave(file.path(OUT, paste0(name,".png")), p, width=w, height=h, dpi=600)
  cat("  :", name, "\n")
}

cat(" …\n")
o <- readRDS("DATA_ROOT")


fp <- function(genes, name, title, ncol=4, w=7.2, h=NULL) {
  g <- intersect(genes, rownames(o))
  if (!length(g)) { cat("  （）:", name, "\n"); return(invisible(NULL)) }
  nr <- ceiling(length(g)/ncol)
  if (is.null(h)) h <- 1.85*nr + 0.55
  pl <- FeaturePlot(o, features=g, reduction="umap", min.cutoff="q5", max.cutoff="q95",
                    pt.size=.5, combine=FALSE)
  pl <- lapply(seq_along(pl), function(i)
    pl[[i]] + scale_color_gradientn(colours=YEL, name=g[i]) + theme_pub +
      theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
            plot.title=element_blank(), legend.position="right",
            legend.key.width=unit(.18,"cm"), legend.key.height=unit(.75,"cm")))
  p <- wrap_plots(pl, ncol=ncol) +
    plot_annotation(title=title, theme=theme(plot.title=element_text(size=10, face="bold", hjust=.5)))
  save_panel(p, name, w, h)
}

PANEL <- list(
  "F10_FeaturePlot_inhibitory_receptors" = list(
    c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","CD160","VSIR","CD244","LAIR1"),
    "Inhibitory receptors"),
  "F11_FeaturePlot_transcription_factors" = list(
    c("TOX","TOX2","NR4A1","NR4A2","NR4A3","BATF","IRF4","PRDM1","ID2","EOMES","TBX21","IKZF2"),
    "Transcription factors"),
  "F12_FeaturePlot_cytotoxic_cytokines" = list(
    c("GZMB","PRF1","IFNG","CCL5","CXCL13"),
    "Cytotoxic effectors and cytokines"),
  "F13_FeaturePlot_surface_markers" = list(
    c("TNFRSF18","TNFRSF9","ITGAE","LAYN","ENTPD1","SLAMF7","CD27","CD28"),
    "Surface and co-signalling markers"),
  "F14_FeaturePlot_memory_progenitor" = list(
    c("TCF7","SELL","CCR7","IL7R","LEF1"),
    "Memory / progenitor markers"),
  "F15_FeaturePlot_cytotoxic_core" = list(
    c("CD8A","CD8B","GZMK","GZMA","NKG7","MKI67","TOP2A","CD69"),
    "Core CD8 and proliferation markers")
)
for (nm in names(PANEL)) {
  a <- PANEL[[nm]]
  fp(a[[1]], nm, a[[2]])
}
cat("\n，:", OUT, "\n")
