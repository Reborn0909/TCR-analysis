
packages <- c("Seurat", "SeuratObject", "harmony", "DoubletFinder", "scRepertoire", "sctransform", "Signac", "Matrix", "edgeR", "limma")
rows <- lapply(packages, function(pkg) data.frame(software = pkg,
    version = tryCatch(as.character(packageVersion(pkg)), error = function(e) "NOT_INSTALLED"),
    scope = "Current installed audit environment; historical analysis version NOT established"))
rows <- c(list(data.frame(software = "R", version = as.character(getRversion()), scope = "Current audit environment")), rows)
write.csv(do.call(rbind, rows), "results/qc_current_software_versions.csv", row.names = FALSE)
print(do.call(rbind, rows))
