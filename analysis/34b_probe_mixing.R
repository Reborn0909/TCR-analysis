
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
setwd("DATA_ROOT")
o <- readRDS("DATA_ROOT")
o$donor <- unname(ifelse(grepl("-1$", as.character(o$sample)), "Donor 1", "Donor 2"))
cat("  :", ncol(o), "· :", paste(sprintf("%s %.3f", names(table(o$donor)), as.numeric(table(o$donor))/ncol(o)), collapse=" | "), "\n")
for (gn in names(o@graphs)) {
  g <- o@graphs[[gn]]
  rs <- Matrix::rowSums(g > 0)
  cat(sprintf("   %-12s  min=%d =%.0f max=%d\n", gn, min(rs), median(rs), max(rs)))
}
g <- o@graphs[["RNA_nn"]]
M <- as(g, "dgCMatrix"); M@x <- rep(1, length(M@x)); M <- M / Matrix::rowSums(M)
D <- as.matrix(Matrix(model.matrix(~ 0 + factor(o$donor)), sparse=TRUE))
frac <- as.matrix(M %*% D); colnames(frac) <- levels(factor(o$donor))
own <- ifelse(o$donor == "Donor 1", frac[,1], frac[,2])
other <- 1 - own
p <- as.numeric(table(o$donor))/ncol(o); exp_same <- sum(p^2); exp_other <- 1 - exp_same
cat(sprintf("\n   p = %.3f / %.3f →   %.4f ·  %.4f\n", p[1], p[2], exp_same, exp_other))
cat("\n  ：\n")
cat(sprintf("    A               = %.4f\n", mean(other)))
cat(sprintf("    B (A - )/(1 - ) = %.4f\n", (mean(other) - exp_other)/(1 - exp_other)))
nd <- rowSums(frac > 0)
cat(sprintf("    C  (k=20)            = %.4f\n", mean(nd)))
cat(sprintf("    D (C-1)/(2-1)                      = %.4f\n", (mean(nd)-1)/1))
cat(sprintf("    E                = %.4f\n", mean(own)))
cat(sprintf("    F 1 - E                            = %.4f\n", 1 - mean(own)))
cat("\n   donor  0.67 —— ：\n")
cands <- c(A=mean(other), B=(mean(other)-exp_other)/(1-exp_other), C=mean(nd), D=(mean(nd)-1), E=mean(own), F=1-mean(own))
print(round(sort(abs(cands - 0.67)), 4))
