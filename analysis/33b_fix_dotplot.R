
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(scales)})
setwd(Sys.getenv("DATA_ROOT", "."))
OUT <- "figures/teacher_panels_umap_dot"
PANELS <- list(
 "Exhaustion" = c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT","TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1"),
 "Activation / effector function" = c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1","KLRD1","KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2","ID2","ID3","PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1","NR4A2","NR4A3","PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB","JUND","STAT1","STAT2","STAT5A","STAT6","STAT4","EOMES"),
 "Cytotoxicity" = c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1","IFNG","TNF","SERPINB1","SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7","CSTB","LAMP1","LAMP3","CAPN2"))
GRP <- c("H-Trp","E-Trp","E-DP")
RAMP <- c("#EFEAE0", "#B5C1B4", "#8C5A5A")
tag <- c("Exhaustion"="exhaustion","Activation / effector function"="activation_effector","Cytotoxicity"="cytotoxicity")

o <- readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
o$group <- factor(o$group, levels=GRP)
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
o$dc <- factor(paste(o$donor, as.character(o$group), sep=" / "),
               levels=as.vector(t(outer(c("Donor 1","Donor 2"), GRP, paste, sep=" / "))))


allg <- unique(unlist(PANELS))
mm <- DotPlot(o, features=allg, group.by="dc", scale=TRUE)$data
lim <- max(abs(range(mm$avg.exp.scaled, na.rm=TRUE)))
lim <- ceiling(lim * 10) / 10
cat(sprintf("  z  ±%.2f →  ±%.1f\n", lim, lim))

for (pn in names(PANELS)) {
  g <- intersect(PANELS[[pn]], rownames(o))
  p <- DotPlot(o, features=g, group.by="dc", cols=RAMP, scale=TRUE, dot.scale=4.6) +
    scale_colour_gradientn(colours=RAMP, limits=c(-lim, lim), oob=scales::squish,
                           name="Average expression\n(scaled)") +
    guides(size=guide_legend(title="Percent expressed", order=2),
           colour=guide_colourbar(order=1, barheight=unit(28,"mm"))) +
    RotatedAxis() +
    labs(title=sprintf("%s: expression across donor and condition", pn),
         subtitle="dot size = fraction of cells expressing; colour = scaled mean expression",
         x=NULL, y=NULL) +
    theme_bw(base_size=9.5) +
    theme(panel.grid.minor=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=7.5),
          plot.title=element_text(face="bold", size=11),
          plot.subtitle=element_text(size=8.5, margin=margin(b=6)),
          legend.position="right", legend.box="vertical",
          plot.margin=margin(6, 8, 6, 6))
  f <- sprintf("D%d_DotPlot_%s.pdf", which(names(PANELS)==pn), tag[[pn]])
  ggsave(file.path(OUT, f), p, width=max(6.6, .26*length(g)+2.6), height=3.4, device=cairo_pdf)
  cat("  :", f, "\n")
}
cat("  \n")
