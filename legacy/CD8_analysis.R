




setwd(Sys.getenv("DATA_ROOT", "."))


dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)


packages_needed <- c("Seurat", "Signac", "tidyverse", "patchwork",
                     "ggplot2", "scRepertoire", "future", "data.table")

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


plan("multisession", workers = 8)
options(future.globals.maxSize = 8000 * 1024^2)


set.seed(42)

cat("，：", getwd(), "\n")






data_dir <- file.path(getwd(), "data")
setwd(data_dir)


zenodo_base <- "https://zenodo.org/records/8240488/files"


files_to_download <- c(
  "PBMC_vaccine_CITE.rds",
  "PBMC_vaccine_ASAP.rds",
  "PBMC_vaccine_ECCITE_TCR.rds",
  "antigen_module_genes.rds",
  "antigen_module_peaks.rds",
  "PBMC_vaccine_ASAP_fragments.tsv.gz",
  "PBMC_vaccine_ASAP_fragments.tsv.gz.tbi"
)


download_with_resume <- function(url, destfile) {
  if (file.exists(destfile)) {
    cat("，：", destfile, "\n")
  } else {
    cat("：", destfile, "\n")
    tryCatch({
      download.file(url, destfile, mode = "wb", method = "auto")
    }, error = function(e) {
      cat("，：", url, "\n")
    })
  }
}

for (file in files_to_download) {
  url <- paste(zenodo_base, file, sep = "/")
  download_with_resume(url, file)
}

cat("\n")


setwd(Sys.getenv("DATA_ROOT", "."))