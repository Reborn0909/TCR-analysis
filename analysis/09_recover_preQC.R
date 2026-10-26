









suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "DATA_ROOT"
SAMP <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("E-CD8DP-1", "E-CD8DP-2")
OBJ  <- "DATA_ROOT"

cat("==  ==\n"); cat(R.version.string, "|", .Platform$OS.type, "\n"); cat(":", .libPaths()[1], "\n")
cat(":", ROOT, "\n:", paste(SAMP, collapse = ", "), "\n\n")
if (!dir.exists(ROOT)) stop("：", ROOT)

read_one <- function(s) {
  cands <- c(file.path(ROOT, s, "filtered_feature_bc_matrix"),
             file.path(ROOT, s),
             file.path(ROOT, paste0("run_count_", s, "-10X5"), "filtered_feature_bc_matrix"),
             file.path(ROOT, paste0(s, "-10X5"), "filtered_feature_bc_matrix"),
             file.path(ROOT, s, paste0(s, "-10X5"), "filtered_feature_bc_matrix"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) return(NULL)
  bc <- file.path(d, "barcodes.tsv.gz"); if (!file.exists(bc)) bc <- file.path(d, "barcodes.tsv")
  ft <- file.path(d, "features.tsv.gz"); if (!file.exists(ft)) ft <- file.path(d, "features.tsv")
  mx <- file.path(d, "matrix.mtx.gz");   if (!file.exists(mx)) mx <- file.path(d, "matrix.mtx")
  if (!all(file.exists(bc, ft, mx))) return(NULL)
  barcodes <- readLines(bc)
  m <- readMM(mx)
  if (nrow(m) != length(barcodes)) {
    if (ncol(m) == length(barcodes)) m <- t(m) else stop("：", d)
  }
  list(dir = d, n = length(barcodes),
       n_ge200 = sum(colSums(m > 0) >= 200),
       barcodes = barcodes)
}

res <- list()
for (s in SAMP) {
  r <- read_one(s)
  if (is.null(r)) { cat(sprintf("%-14s \n", s)); next }
  cat(sprintf("%-14s barcodes=%d   min.features>=200=%d  =%d\n",
              s, r$n, r$n_ge200, r$n - r$n_ge200))
  cat(sprintf("                %s\n", r$dir))
  res[[s]] <- r
}
if (!length(res)) stop("；。")

cat("\n==  newQC （）==\n")
ok_all <- TRUE

CELLCSV <- "results/main_metadata_verified_from_original.csv"
cell_names_by_sample <- NULL
if (file.exists(CELLCSV) && Sys.getenv("FORCE_OBJ") != "1") {
  md <- read.csv(CELLCSV, row.names = 1, check.names = FALSE)
  if (all(c("sample") %in% names(md))) {
    cell_names_by_sample <- split(rownames(md), md$sample)
    cat("：", CELLCSV, "（", nrow(md), " ）\n", sep = "")
  }
}
if (!is.null(cell_names_by_sample)) {
  for (s in names(res)) {
    post <- cell_names_by_sample[[s]]
    if (is.null(post)) post <- character(0)
    if (!length(post)) { cat(sprintf("%-14s \n", s)); next }


    bc0 <- res[[s]]$barcodes
    p0  <- sub(paste0("^", s, "_"), "", post)
    inside <- sum(p0 %in% bc0)
    if (inside == 0 && length(post) > 0) {
      bc1 <- sub("-1$", "", bc0); p1 <- sub("-1$", "", p0)
      inside <- sum(p1 %in% bc1)
      if (inside > 0) cat(sprintf("%-14s ： -1 \n", s))
    }
    good <- inside == length(post)
    ok_all <- ok_all && good
    cat(sprintf("%-14s QC %d   %d   %d  %s\n", s, length(post), inside,
                length(post) - inside, if (good) " OK" else "，"))
  }
} else if (file.exists(OBJ) && Sys.getenv("SKIP_OBJ") != "1") {
  suppressPackageStartupMessages(library(Seurat))
  o <- readRDS(OBJ)
  for (s in names(res)) {
    post <- colnames(o)[o$sample == s]
    if (!length(post)) { cat(sprintf("%-14s \n", s)); next }
    bc0 <- res[[s]]$barcodes
    p0  <- sub(paste0("^", s, "_"), "", post)
    inside <- sum(p0 %in% bc0)
    good <- inside == length(post)
    ok_all <- ok_all && good
    cat(sprintf("%-14s QC %d   %d   %d  %s\n", s, length(post), inside,
                length(post) - inside, if (good) " OK" else "，"))
  }
} else cat("（）\n")

out <- data.frame(
  sample = names(res),
  pre_QC_cells = vapply(res, function(x) x$n, 0L),
  passing_min_features_200 = vapply(res, function(x) x$n_ge200, 0L),
  removed_by_min_features_200 = vapply(res, function(x) x$n - x$n_ge200, 0L),
  matrix_dir = vapply(res, function(x) x$dir, ""),
  stringsAsFactors = FALSE)
f <- "results/main_preQC_cell_counts_recovered.csv"
write.csv(out, f, row.names = FALSE)
cat("\n：", f, "\n"); print(out)
if (!ok_all) cat("\n：，，。\n")
