












suppressPackageStartupMessages({library(Matrix); library(Seurat); library(FNN)})
ROOT <- "F:/Jpan_20261007/Trp_10X5/Trp_10X5"
SAMPLES <- c("E-CD8DP-1","E-CD8DP-2","E-CD8Trp-1","E-CD8Trp-2","H-CD8Trp-1","H-CD8Trp-2")
set.seed(20261008)
N_SIM <- 1500; NPC <- 30; K <- 30
exp_rate <- function(n) approx(c(500,1000,2000,5000,10000,15000),
                               c(0.4,0.8,1.6,3.9,7.6,11.0), xout = n, rule = 2)$y / 100
res <- list()
for (s in SAMPLES) {
  d <- file.path(ROOT, sprintf("run_count_%s-10X5", s), "filtered_feature_bc_matrix")
  if (!dir.exists(d)) { cat(sprintf("%-12s \n", s)); next }
  m <- Read10X(d)
  feat <- colSums(m > 0); cnt <- colSums(m)
  mtg <- grep("^MT-", rownames(m), value = TRUE)
  mt <- if (length(mtg)) colSums(m[mtg, , drop = FALSE]) / cnt * 100 else rep(0, ncol(m))
  hbg <- grep("^HB[AB]", rownames(m), value = TRUE)
  hb <- if (length(hbg)) colSums(m[hbg, , drop = FALSE]) / cnt * 100 else rep(0, ncol(m))
  k <- feat > 200 & feat < 10000 & cnt < 80000 & mt < 10 & hb < 0.35
  X <- m[, k, drop = FALSE]; X <- X[Matrix::rowSums(X > 0) >= 3, , drop = FALSE]
  n <- ncol(X)
  i1 <- sample.int(n, N_SIM, replace = TRUE); i2 <- sample.int(n, N_SIM, replace = TRUE)
  Sim <- X[, i1, drop = FALSE] + X[, i2, drop = FALSE]
  colnames(Sim) <- paste0("sim", seq_len(N_SIM))
  obj <- CreateSeuratObject(counts = cbind(X, Sim)); rm(X, Sim, m); gc(verbose = FALSE)
  obj <- NormalizeData(obj, verbose = FALSE); obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE); obj <- RunPCA(obj, npcs = NPC, verbose = FALSE)
  emb <- Embeddings(obj, "pca")[, 1:NPC, drop = FALSE]
  nn <- get.knnx(emb, emb[seq_len(n), , drop = FALSE], k = K)$nn.index
  score <- rowMeans(matrix(nn > n, nrow = n))
  er <- exp_rate(n)
  res[[s]] <- data.frame(
    sample = s, analyzed_cells = n,
    expected_rate_pct_10x = round(er * 100, 2),
    score_median = round(median(score), 4),
    score_p90 = round(quantile(score, .90), 4),
    score_p99 = round(quantile(score, .99), 4),
    frac_ge_0.20 = round(mean(score >= 0.20) * 100, 2),
    frac_ge_0.33 = round(mean(score >= 0.33) * 100, 2),
    frac_ge_0.50 = round(mean(score >= 0.50) * 100, 2))
  cat(sprintf("%-12s n=%-6d %.2f%%  =%.3f p90=%.3f  ≥0.20:%.1f%% ≥0.33:%.1f%% ≥0.50:%.1f%%\n",
      s, n, er*100, median(score), quantile(score,.90),
      mean(score>=.20)*100, mean(score>=.33)*100, mean(score>=.50)*100))
  rm(obj, emb); gc(verbose = FALSE)
}
out <- do.call(rbind, res)
write.csv(out, "results/main_doublet_assessment.csv", row.names = FALSE)
cat("\n results/main_doublet_assessment.csv（）\n"); print(out)
