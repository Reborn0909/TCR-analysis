

suppressPackageStartupMessages({ library(SeuratObject); library(data.table) })
set.seed(20261004)

OBJ <- Sys.getenv("OBJECT", file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))
OUT <- Sys.getenv("SCORES_DIR", "data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console()

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

msg("[1/5]  RDS（2.6 GB，）...")
obj <- readRDS(OBJ)
msg("      class:", paste(class(obj), collapse="/"), "· cells:", ncol(obj), "· genes:", nrow(obj))


DefaultAssay(obj) <- "RNA"
dat <- tryCatch(GetAssayData(obj, assay="RNA", layer="data"),
                error=function(e) GetAssayData(obj, assay="RNA", slot="data"))
msg("      data layer:", class(dat)[1], nrow(dat), "x", ncol(dat))


avg <- Matrix::rowMeans(dat)
msg("[2/5]  24 bins（Seurat ）")
nbin <- 24; ctrl <- 100
grp <- cut(avg, breaks = nbin, labels = FALSE)
names(grp) <- rownames(dat)
tab <- table(grp)
msg("       bin :", paste(head(as.integer(tab), 6), collapse=","), "...")

control_sets <- list()
for (nm in names(PANELS)) {
  fs <- intersect(PANELS[[nm]], rownames(dat))
  msg("", nm, ": ", length(fs), "/", length(PANELS[[nm]]), "")

  cg <- character(0)
  for (b in unique(grp[fs])) {
    pool <- setdiff(names(grp)[grp == b], fs)
    if (length(pool) == 0) next
    cg <- c(cg, sample(pool, min(ctrl, length(pool))))
  }
  control_sets[[nm]] <- unique(cg)
  msg("        :", length(control_sets[[nm]]))
}

msg("[3/5] （ - ）")
out <- data.table(cell = colnames(dat))
for (nm in names(PANELS)) {
  fs <- intersect(PANELS[[nm]], rownames(dat))
  cg <- intersect(control_sets[[nm]], rownames(dat))
  sc <- Matrix::colMeans(dat[fs, , drop=FALSE]) - Matrix::colMeans(dat[cg, , drop=FALSE])
  out[[nm]] <- as.numeric(sc)
  msg("", nm, "  score: mean=", round(mean(sc),4), " sd=", round(sd(sc),3))
}

msg("[4/5] ")
md <- obj@meta.data
msg("      meta :", paste(head(colnames(md), 12), collapse=", "))
smp_col <- intersect(c("sample","orig.ident","Sample"), colnames(md))[1]
smp <- if (!is.na(smp_col)) as.character(md[[smp_col]]) else colnames(obj)
out[, sample := smp]
out[, group := ifelse(grepl("^H-", smp), "SS",
               ifelse(grepl("^E-CD8Trp", smp), "Pre-sel",
               ifelse(grepl("^E-CD8DP", smp), "Post-sel", NA)))]
out[, donor := ifelse(grepl("-1$", smp), "Donor1",
               ifelse(grepl("-2$", smp), "Donor2", NA))]
msg("[5/5] ")
print(table(out$group, out$donor, useNA="ifany"))
fwrite(out, file.path(OUT, "per_cell_scores.csv.gz"))
msg(" → ", file.path(OUT, "per_cell_scores.csv.gz"), " (", nrow(out), " )")
