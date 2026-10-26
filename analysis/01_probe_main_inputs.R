args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('Usage: Rscript 01_probe_main_inputs.R seurat.rds combined_tcr.rds output_dir')
suppressPackageStartupMessages(library(Seurat))
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
cat('Run:', format(Sys.time(), tz='Asia/Shanghai'), '\n')
cat('Input:', args[1], '\n')
tcr <- readRDS(args[2])
cat('Combined TCR elements:\n'); print(names(tcr))
for (n in names(tcr)) {
 cat('\nTCR:', n, 'n=', nrow(tcr[[n]]), '\n')
 print(names(tcr[[n]])); print(utils::head(tcr[[n]], 1))
}
obj <- readRDS(args[1])
cat('\nSeurat version stored:', as.character(obj@version), '\n')
cat('Assays/layers:\n'); for (a in Assays(obj)) {print(a); print(Layers(obj[[a]])); for (l in Layers(obj[[a]])) cat(l, paste(dim(LayerData(obj, assay=a, layer=l)), collapse=' x '), '\n')}
cat('\nCommands:\n')
for (n in names(obj@commands)) {cat('\n',n,'\n'); print(obj@commands[[n]]@params)}
md <- obj[[]]
write.csv(md, file.path(args[3], 'main_metadata_verified_from_original.csv'), row.names=TRUE)
qc <- do.call(rbind, lapply(split(md,md$sample),function(d) data.frame(sample=d$sample[1],group=d$group[1],n_post=nrow(d),median_genes=median(d$nFeature_RNA),median_UMI=median(d$nCount_RNA),median_pct_mt=median(d$percent.mt),CTstrict_nonmissing=sum(!is.na(d$CTstrict)),CTstrict_missing_as_single=sum(is.na(d$CTstrict)&d$cloneSize_real==1))))
write.csv(qc, file.path(args[3], 'main_qc_verified_from_original.csv'),row.names=FALSE)
cat('\nOriginal QC: \n');print(qc)
cat('\nSession:\n');print(sessionInfo())
