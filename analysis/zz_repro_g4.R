

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
B <- "DATA_ROOT"
long <- read.csv(paste0(B, "expr.csv"), fileEncoding="UTF-8-BOM")
agg  <- read.csv(paste0(B, "delta.csv"), fileEncoding="UTF-8-BOM")
GRP  <- c("H-Trp","E-Trp","E-DP")
HEAT <- c("#D4CFC5","#B5C1B4","#8C5A5A")
hm <- long %>% mutate(col = paste(donor, group, sep=" · ")) %>% select(gene, col, mean_log1p)
ord <- agg %>% arrange(panel, desc(mean_delta)) %>% pull(gene)
hm$gene <- factor(hm$gene, levels=rev(ord))
hm$col  <- factor(hm$col, levels=as.vector(t(outer(c("Donor 1","Donor 2"), GRP, paste, sep=" · "))))
hm <- hm %>% group_by(gene) %>% mutate(z = as.numeric(scale(mean_log1p))) %>% ungroup()
cat("", nrow(hm), "· z ", sum(!is.finite(hm$z)), "\n")
cat("  col  NA :", sum(is.na(hm$col)), "· gene  NA :", sum(is.na(hm$gene)), "\n")
cat("  z :", sprintf("%.3f ~ %.3f", min(hm$z, na.rm=TRUE), max(hm$z, na.rm=TRUE)), "\n")
p <- ggplot(hm, aes(col, gene, fill=z)) + geom_tile(colour="white", linewidth=.25) +
  scale_fill_gradientn(colours=HEAT, limits=c(-2.4,2.4), oob=scales::squish, name="row z") +
  labs(title="repro", x=NULL, y=NULL) + theme_minimal(base_size=6)
ggsave("DATA_ROOT",
       p, width=4.6, height=11.5, device=cairo_pdf)
cat("   repro.pdf\n")
