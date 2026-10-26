




suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
RES <- "DATA_ROOT"
theme_pub <- theme_minimal(base_size = 8) +
  theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 7.5, hjust = 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7, color="black"),
        panel.border = element_rect(linewidth=.3, fill=NA, color="black"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.position = "bottom")
theme_set(theme_pub)

CONTRASTS <- list(
  "E_Trp_vs_H_Trp" = c("E-Trp vs H-Trp", "E-Trp minus H-Trp"),
  "E_DP_vs_E_Trp"  = c("E-DP vs E-Trp",  "E-DP minus E-Trp"),
  "E_DP_vs_H_Trp"  = c("E-DP vs H-Trp",  "E-DP minus H-Trp"))

for (nm in names(CONTRASTS)) {
  f <- file.path(RES, sprintf("main_pseudobulk_DE_%s.csv", nm))
  if (!file.exists(f)) { cat("  :", nm, "\n"); next }
  z <- read_csv(f, show_col_types = FALSE) %>%
    mutate(sig = ifelse(FDR < 0.05 & abs(logFC) > 1, "FDR < 0.05 and |log2FC| > 1", "Other"),
           neglogFDR = -log10(pmax(FDR, 1e-300)))
  z$sig <- factor(z$sig, levels = c("FDR < 0.05 and |log2FC| > 1", "Other"))
  p <- ggplot(z, aes(logFC, neglogFDR, colour = sig)) +
    geom_point(size = .5, alpha = .60) +
    scale_colour_manual(values = c("FDR < 0.05 and |log2FC| > 1" = "#B05A5A", "Other" = "#B8B8B8")) +
    labs(title = sprintf("Donor-paired pseudobulk: %s", CONTRASTS[[nm]][1]),
         subtitle = sprintf("%s   |   edgeR TMM + quasi-likelihood, ~ donor + group; two donors (exploratory)",
                            CONTRASTS[[nm]][2]),
         x = expression(log[2]~"fold change"), y = expression(-log[10]~"FDR"))
  out <- file.path(RES, sprintf("main_pseudobulk_volcano_%s", nm))
  cairo_pdf(paste0(out, ".pdf"), width = 6.6, height = 5.0); print(p); dev.off()
  ggsave(paste0(out, ".png"), p, width = 6.6, height = 5.0, dpi = 600)
  cat("  :", nm, " n=", nrow(z), " =", sum(z$sig == "FDR < 0.05 and |log2FC| > 1"), "\n", sep="")
}
