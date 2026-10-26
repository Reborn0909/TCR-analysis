





suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggsignif); library(patchwork)
})

OUT <- Sys.getenv("GHD_OUT", "/mnt/f// for GHD/main_figure")
IN  <- file.path(OUT, "per_cell_scores.csv.gz")
if (!file.exists(IN)) stop(" per_cell_scores.csv.gz， 01_extract_scores.R")

dt <- fread(IN)
dt[, group := factor(group, levels = c("SS","Pre-sel","Post-sel"))]
dt[, donor := factor(donor, levels = c("Donor 1","Donor 2"))]

PAL <- c("SS" = "#8FA9B8", "Pre-sel" = "#C08A8A", "Post-sel" = "#D0A97E")
PANEL_TITLES <- c(Exhaustion = "Exhaustion", Activation_effector_function = "Activation / Effector function",
                  Cytotoxicity = "Cytotoxicity")

make_panel <- function(pn, title_txt) {
  sub <- dt[!is.na(get(pn)) & !is.na(group)]
  rng <- range(sub[[pn]], na.rm = TRUE); pad <- diff(rng) * 0.30
  ymax <- max(sub[[pn]], na.rm = TRUE) + pad

  cmp <- list(c("SS","Pre-sel"), c("Pre-sel","Post-sel"), c("SS","Post-sel"))
  pos <- c(ymax, ymax + diff(rng)*0.07, ymax + diff(rng)*0.14)
  ggplot(sub, aes(x = group, y = get(pn), fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.85, colour = "grey30", linewidth = 0.35) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9,
                 colour = "grey20", linewidth = 0.35) +
    ggsignif::geom_signif(comparisons = cmp, annotations = c("****","****","****"),
                          y_position = pos, tip_length = 0.012,
                          vjust = 0.35, textsize = 3.2, linewidth = 0.35) +
    scale_fill_manual(values = PAL) +
    labs(title = title_txt, x = NULL, y = "Score") +
    coord_cartesian(ylim = c(rng[1] - pad, max(pos) + pad*0.4)) +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 9, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          axis.title.y = element_text(size = 9),
          plot.margin = margin(6, 8, 4, 4))
}

ps <- list(
  make_panel("Exhaustion", "Exhaustion"),
  make_panel("Activation_effector_function", "Activation / Effector function"),
  make_panel("Cytotoxicity", "Cytotoxicity")
)
combined <- (ps[[1]] | ps[[2]] | ps[[3]]) + plot_layout(ncol = 3)

ggsave(file.path(OUT, "main_figure_module_scores.pdf"), combined,
       width = 7.4, height = 2.9, device = cairo_pdf)
ggsave(file.path(OUT, "main_figure_module_scores.png"), combined,
       width = 7.4, height = 2.9, dpi = 600)
message(" → ", file.path(OUT, "main_figure_module_scores.pdf"))
