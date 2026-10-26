## Depth-matched TCR diversity: does the SS -> Pre-sel -> Post-sel increase
## survive when every sample is held to the same per-cell VDJ depth?
ITER <- 2000; SEED <- 20261004
D <- "/tmp/ghd_productive/downstream_clean"
out <- "/tmp/ghd_productive/depth_matched"
dir.create(out, showWarnings = FALSE)
shannon <- function(k) { p <- table(k)/length(k); -sum(p*log(p)) }
ORDER <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")
GRP <- c("H-CD8Trp-1"="SS","E-CD8Trp-1"="Pre","E-CD8DP-1"="Post",
         "H-CD8Trp-2"="SS","E-CD8Trp-2"="Pre","E-CD8DP-2"="Post")
DON <- c("H-CD8Trp-1"="Donor1","E-CD8Trp-1"="Donor1","E-CD8DP-1"="Donor1",
         "H-CD8Trp-2"="Donor2","E-CD8Trp-2"="Donor2","E-CD8DP-2"="Donor2")

for (nm in c("CTgene_like","CTstrict_like")) {
  ck <- readRDS(file.path(D, paste0("cellkeys_", nm, ".rds")))
  cat("\n=========================", nm, "=========================\n")
  rows <- list()
  for (fl in c(1, 5, 10, 15, 20)) {
    sub <- ck[ck$umis >= fl, ]
    cnt <- table(factor(sub$sample, levels = ORDER))
    n <- min(cnt)
    if (n < 200) { cat("floor", fl, ": min cells", n, "-> skipped\n"); next }
    set.seed(SEED)
    for (s in ORDER) {
      g <- sub[sub$sample == s, ]
      v <- numeric(ITER)
      for (i in seq_len(ITER)) v[i] <- shannon(g$key[sample(seq_len(nrow(g)), n, replace = FALSE)])
      rows[[length(rows)+1]] <- data.frame(
        floor_umis = fl, sample = s, group = GRP[s], donor = DON[s],
        eligible = nrow(g), n_downsample = n,
        raw_Shannon = shannon(g$key), downsample_median = median(v),
        lower = quantile(v,.025,names=FALSE), upper = quantile(v,.975,names=FALSE))
    }
    cat(sprintf("floor >= %2d UMIs: n_downsample = %d\n", fl, n))
  }
  R <- do.call(rbind, rows)
  write.csv(R, file.path(out, paste0("depth_matched_", nm, ".csv")), row.names = FALSE)
  cat("\n--- median downsampled Shannon by depth floor ---\n")
  m <- reshape(R[,c("floor_umis","sample","downsample_median")], idvar="floor_umis",
               timevar="sample", direction="wide")
  print(m, row.names = FALSE)
  cat("\n--- donor trend (is Pre>SS and Post>Pre at each floor?) ---\n")
  for (fl in sort(unique(R$floor_umis))) for (d in c("Donor1","Donor2")) {
    q <- R[R$floor_umis==fl & R$donor==d, ]
    ss <- q$downsample_median[q$group=="SS"]; pre <- q$downsample_median[q$group=="Pre"]; po <- q$downsample_median[q$group=="Post"]
    cat(sprintf("  floor %2d %s: SS=%.4f Pre=%.4f Post=%.4f  Pre>SS:%s Post>Pre:%s\n",
                fl, d, ss, pre, po, pre>ss, po>pre))
  }
}
