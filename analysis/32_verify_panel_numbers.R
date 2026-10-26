








suppressPackageStartupMessages({library(Seurat); library(Matrix); library(dplyr); library(tidyr)})
setwd(Sys.getenv("DATA_ROOT", "."))

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

o <- readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
o$group <- factor(o$group, levels=c("H-Trp","E-Trp","E-DP"))
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
cat("  :", ncol(o), " ·  assay:", DefaultAssay(o), "\n")


m <- GetAssayData(o, assay="RNA", layer="data")[ALLG, , drop=FALSE]
res <- list()
for (d in c("Donor 1","Donor 2")) for (g in c("H-Trp","E-Trp","E-DP")) {
  idx <- which(o$donor == d & as.character(o$group) == g)
  res[[length(res)+1]] <- data.frame(donor=d, group=g, n_B=length(idx),
                                     gene=rownames(m), mean_B=Matrix::rowMeans(m[, idx, drop=FALSE]))
}
b <- do.call(rbind, res)


a <- read.csv("results/teacher_panel_donor_expression.csv", fileEncoding="UTF-8-BOM")
mrg <- merge(a, b, by=c("donor","group","gene"))
mrg$diff <- abs(mrg$mean_log1p - mrg$mean_B)
mx <- max(mrg$diff, na.rm=TRUE)
cat(sprintf("   %d （ %d）·  %.3e\n", nrow(mrg), 76*6, mx))
cat(sprintf("  : %s\n", all(mrg$n_cells == mrg$n_B)))



w2 <- mrg %>% select(gene, donor, group, mean_B) %>% tidyr::pivot_wider(names_from=group, values_from=mean_B)
d3 <- w2 %>% mutate(delta = `E-Trp` - `H-Trp`) %>% group_by(gene) %>%
  summarise(d1=delta[donor=="Donor 1"], d2=delta[donor=="Donor 2"], .groups="drop") %>%
  mutate(up=d1>0 & d2>0, down=d1<0 & d2<0)
cat(sprintf("    %d ·  %d ·  %d\n",
            sum(d3$up), sum(d3$down), sum(!d3$up & !d3$down)))
cat(sprintf("   LAG3 Δ: D1 %+0.4f  D2 %+0.4f\n",
            d3$d1[d3$gene=="LAG3"], d3$d2[d3$gene=="LAG3"]))



if (nrow(mrg) != 76*6) {
  cat(sprintf("  ** %d ≠ 456 —— merge  ✗**\n", nrow(mrg)))
  quit(status = 1)
}
if (mx < 1e-9 && all(mrg$n_cells == mrg$n_B)) {
  cat("  ** ✓**（ vs FetchData，）\n")
  write.csv(data.frame(check="matrix rowMeans vs FetchData mean", rows=nrow(mrg),
                       max_abs_diff=mx, up=sum(d3$up), down=sum(d3$down),
                       lag3_d1=d3$d1[d3$gene=="LAG3"], lag3_d2=d3$d2[d3$gene=="LAG3"],
                       verdict="identical"),
            "manifests/teacher_panel_recompute_check.csv", row.names=FALSE)
} else {
  cat("  ** ✗ ", mx, "**\n"); quit(status=1)
}
