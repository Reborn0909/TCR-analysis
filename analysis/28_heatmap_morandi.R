




suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
setwd("DATA_ROOT")
OUT <- "figures/panels_en"
RAMP <- c("#D4CFC5", "#B5C1B4", "#8C5A5A")
SAMP <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")
theme_pub <- theme_minimal(base_size=7) +
  theme(plot.title=element_text(size=8.5,face="bold",hjust=.5,margin=margin(b=3)),
        plot.subtitle=element_text(size=6.8,hjust=.5),
        axis.text.x=element_text(size=6,angle=45,hjust=1,color="black"),
        axis.text.y=element_text(size=6,color="black"),
        legend.title=element_text(size=7), legend.text=element_text(size=6),
        legend.key.width=unit(.28,"cm"), legend.key.height=unit(1.5,"cm"),
        panel.grid=element_blank(),
        panel.border=element_rect(linewidth=.3,fill=NA,color="black"))
theme_set(theme_pub)
draw <- function(csv, name, title, w, h) {
  d <- read_csv(csv, show_col_types=FALSE)
  d$sample <- factor(d$sample, levels=SAMP)

  ord <- unique(d$gene)
  d$gene <- factor(d$gene, levels=rev(ord))
  p <- ggplot(d, aes(sample, gene, fill=z)) +
    geom_tile(colour="white", linewidth=.35) +
    scale_fill_gradientn(colours=RAMP, name="Row z-score", limits=c(-2,2), oob=scales::squish) +
    labs(title=title, subtitle="TMM-normalised logCPM; each gene standardised across the six pseudobulk libraries",
         x=NULL, y=NULL)
  cairo_pdf(file.path(OUT, paste0(name,".pdf")), width=w, height=h); print(p); dev.off()
  ggsave(file.path(OUT, paste0(name,".png")), p, width=w, height=h, dpi=600)
  cat("  :", name, " :", length(ord), "\n")
}
draw("results/main_predefined_key_genes_pseudobulk_heatmap_source.csv",
     "P11_predefined_key_genes_heatmap_morandi", "Predefined genes: donor-level expression display", 4.0, 5.4)
draw("results/main_G_pseudobulk_selected_candidates_heatmap_source.csv",
     "P12_candidate_genes_heatmap_morandi", "G candidate replacement: selected from paired pseudobulk DE", 4.0, 5.6)
