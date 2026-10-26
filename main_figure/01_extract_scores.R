





suppressPackageStartupMessages({ library(Seurat); library(data.table) })

OBJ <- Sys.getenv("GHD_OBJECT", "/mnt/f// for GHD/seurat_filtered_newQC.rds")
OUT <- Sys.getenv("GHD_OUT", "/mnt/f// for GHD/main_figure")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PANELS <- list(
  Exhaustion = c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT",
                 "TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1"),
  Activation_effector_function = c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1","KLRD1",
                 "KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2","ID2","ID3",
                 "PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1","NR4A2","NR4A3",
                 "PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB","JUND","STAT1","STAT2",
                 "STAT5A","STAT6","STAT4","EOMES"),
  Cytotoxicity = c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1","IFNG","TNF","SERPINB1",
                   "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3",
                   "CST7","CSTB","LAMP1","LAMP3","CAPN2"))

message("[1/4] : ", OBJ)
obj <- readRDS(OBJ)
DefaultAssay(obj) <- "RNA"
message("     cells = ", ncol(obj), " · genes = ", nrow(obj))


meta <- obj@meta.data
if (!"group" %in% colnames(meta)) {
  smp <- if ("sample" %in% colnames(meta)) as.character(meta$sample) else rownames(meta)
  grp <- ifelse(grepl("^H-", smp), "SS",
         ifelse(grepl("^E-CD8Trp", smp), "Pre-sel",
         ifelse(grepl("^E-CD8DP", smp), "Post-sel", NA)))
  obj$group <- factor(grp, levels = c("SS","Pre-sel","Post-sel"))
  obj$sample_short <- smp
}
obj$donor_f <- factor(ifelse(grepl("-1$", as.character(obj$sample_short)), "Donor 1",
                      ifelse(grepl("-2$", as.character(obj$sample_short)), "Donor 2", "NA")),
                      levels = c("Donor 1","Donor 2"))

message("[2/4]  AddModuleScore")
for (nm in names(PANELS)) {
  g <- intersect(PANELS[[nm]], rownames(obj))
  message("", nm, ": ", length(g), "/", length(PANELS[[nm]]), "")
  obj <- AddModuleScore(obj, features = list(g), name = paste0("MS_", nm), seed = 42)
}
score_cols <- paste0("MS_", names(PANELS), "1")
names(score_cols) <- names(PANELS)

message("[3/4] ")
per_cell <- data.table(
  cell = colnames(obj),
  sample = as.character(obj$sample_short),
  group = as.character(obj$group),
  donor = as.character(obj$donor_f)
)
for (nm in names(PANELS)) per_cell[[nm]] <- obj[[ score_cols[[nm]] ]][, 1]
fwrite(per_cell, file.path(OUT, "per_cell_scores.csv.gz"))

message("[4/4] ")
long <- melt(per_cell, id.vars = c("cell","sample","group","donor"),
             variable.name = "panel", value.name = "score")
summ <- long[, .(n_cells = .N,
                 mean_score = mean(score, na.rm = TRUE),
                 median_score = median(score, na.rm = TRUE),
                 q25 = quantile(score, .25, na.rm = TRUE),
                 q75 = quantile(score, .75, na.rm = TRUE)),
             by = .(donor, group, panel)]
fwrite(summ, file.path(OUT, "module_score_donor_summary.csv"))
message(" → ", OUT)
