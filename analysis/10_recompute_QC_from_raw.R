





suppressPackageStartupMessages({library(Matrix); library(Seurat)})
ROOT <- "F:/GHD_Jpan_20261007/Trp_10X5/Trp_10X5"
OBJ  <- "DATA_ROOT"
META <- "DATA_ROOT"
SAMPLES <- c("E-CD8DP-1","E-CD8DP-2","E-CD8Trp-1","E-CD8Trp-2","H-CD8Trp-1","H-CD8Trp-2")
GROUPS  <- c("E-DP","E-DP","E-Trp","E-Trp","H-Trp","H-Trp")
cat("==  ==\n", R.version.string, "|", .Platform$OS.type, "\n:", .libPaths()[1], "\n\n")

qc <- list(); n <- 0
for (i in seq_along(SAMPLES)) {
  s <- SAMPLES[i]
  d <- file.path(ROOT, sprintf("run_count_%s-10X5", s), "filtered_feature_bc_matrix")
  if (!dir.exists(d)) { cat(sprintf("%-12s \n", s)); next }
  m <- Read10X(d)
  n0 <- ncol(m)
  feat <- colSums(m > 0); cnt <- colSums(m);
  mtg <- grep("^MT-", rownames(m), value = TRUE); hbg <- grep("^HB[AB]", rownames(m), value = TRUE)
  mt <- if (length(mtg)) colSums(m[mtg, , drop = FALSE]) / cnt * 100 else rep(0, n0)
  hb <- if (length(hbg)) colSums(m[hbg, , drop = FALSE]) / cnt * 100 else rep(0, n0)
  k <- feat > 200 & feat < 10000 & cnt < 80000 & mt < 10 & hb < 0.35
  cat(sprintf("%-12s barcodes=%-6d =%-6d (%.1f%%)  : feat%d feat%d UMI%d mt%d hb%d\n",
      s, n0, sum(k), 100*sum(k)/n0,
      sum(feat <= 200), sum(feat >= 10000), sum(cnt >= 80000), sum(mt >= 10), sum(hb >= 0.35)))
  qc[[s]] <- data.frame(sample = s, group = GROUPS[i], pre_QC = n0, post_QC = sum(k))
  n <- n + sum(k)
}
cat("\n QC :", n, "\n")


if (file.exists(META)) {
  md <- read.csv(META, row.names = 1, check.names = FALSE)
  cat("\n==  ==\n")
  for (s in names(qc)) {
    obj_n <- sum(md$sample == s)
    cat(sprintf("%-12s  %-6d  %-6d %s\n", s, qc[[s]]$post_QC, obj_n,
                ifelse(qc[[s]]$post_QC == obj_n, "", "")))
  }
}
out <- do.call(rbind, qc)
write.csv(out, "results/main_raw_QC_recomputation.csv", row.names = FALSE)
cat("\n results/main_raw_QC_recomputation.csv\n"); print(out)
