


suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
source_root <- ".."
arg_index <- 1L
while (arg_index <= length(args)) {
  if (args[arg_index] == "--source-root") {
    if (arg_index == length(args)) stop("--source-root requires a directory path")
    source_root <- args[arg_index + 1L]
    arg_index <- arg_index + 2L
  } else if (args[arg_index] %in% c("--full-vaccine", "--sle-final")) {
    arg_index <- arg_index + 1L
  } else stop("Unknown argument: ", args[arg_index])
}
if (!dir.exists(source_root)) stop("Source directory does not exist: ", source_root)
if (all(c("--full-vaccine", "--sle-final") %in% args)) stop("Choose one mode: --full-vaccine or --sle-final")
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
files <- c(vaccine_CD8 = "../data/cd8_cite_final_ref.rds",
           vaccine_CD4 = "../data/cd4_cells_final_ref.rds",
           vaccine_ECCITE = "../data/cd8_eccite_analyzed.rds",
           SLE_GSE135779_metadata = "../results/SLE_GSE135779/cell_metadata.rds")
files <- setNames(file.path(source_root, sub("^\\.\\./", "", files)), names(files))
if ("--full-vaccine" %in% commandArgs(trailingOnly = TRUE)) {
  files <- c(vaccine_CITE_full = file.path(source_root, "data/PBMC_vaccine_CITE.rds"),
             vaccine_ECCITE_full = file.path(source_root, "data/PBMC_vaccine_ECCITE_TCR.rds"))
}
if ("--sle-final" %in% args) {
  files <- c(SLE_GSE135779_CD8_final = file.path(source_root, "results/SLE_GSE135779/CD8_SLE_final.rds"),
             SLE_GSE135779_CD4_final = file.path(source_root, "results/SLE_GSE135779/CD4_SLE_final.rds"))
}
summary_rows <- list()
for (label in names(files)) {
  cat("\nOBJECT", label, files[[label]], "\n")
  if (!file.exists(files[[label]])) { cat("MISSING\n"); next }
  obj <- readRDS(files[[label]])
  cat("class:", paste(class(obj), collapse = ","), "\n")
  if (inherits(obj, "Seurat")) {
    md <- obj@meta.data
    cat("stored Seurat object version:", as.character(obj@version), "\n")
    cat("cells:", nrow(md), "genes:", nrow(obj), "\n")
    cat("assays:", paste(Assays(obj), collapse = ","), "\n")
    cat("reductions:", paste(Reductions(obj), collapse = ","), "\n")
    cat("COMMANDS\n")
    for (nm in names(obj@commands)) {
      cat("\n", nm, "\n")
      print(obj@commands[[nm]]@params)
    }
    if (length(obj@commands)) {
      command_rows <- lapply(names(obj@commands), function(nm) {
        data.frame(object = label, command = nm,
                   params = paste(capture.output(print(obj@commands[[nm]]@params)), collapse = "\n"))
      })
      write.csv(do.call(rbind, command_rows), file.path("results", paste0("external_", label, "_commands.csv")), row.names = FALSE, fileEncoding = "UTF-8")
    }
    object_note <- if (grepl("_full$", label)) "Published processed full object, after original QC; before-QC counts unavailable" else if (grepl("SLE_GSE135779_.*_final$", label)) "Directly read local GSE135779 final subset; not GSE174188 and before-QC counts unavailable" else "Existing processed subset, not raw input or complete PBMC post-QC count"
    summary_rows[[label]] <- data.frame(label = label, source = files[[label]], cells = nrow(md), genes = nrow(obj), stored_object_version = as.character(obj@version), reductions = paste(Reductions(obj), collapse = ";"), note = object_note)
  } else if (is.data.frame(obj)) {
    md <- obj
    cat("metadata cells:", nrow(md), "\n")
  } else {
    str(obj, max.level = 1)
    rm(obj); gc(); next
  }
  cat("metadata columns:", paste(names(md), collapse = ","), "\n")
  metadata_columns <- data.frame(object = label, column = names(md),
                                class = vapply(md, function(x) paste(class(x), collapse = ";"), character(1)),
                                qc_candidate = grepl("nFeature|nCount|percent|mito|doublet|scrublet", names(md), ignore.case = TRUE))
  write.csv(metadata_columns, file.path("results", paste0("external_", label, "_metadata_columns.csv")), row.names = FALSE, fileEncoding = "UTF-8")
  if (all(c("donor", "timepoint") %in% names(md))) {
    md$donor_timepoint <- paste(md$donor, md$timepoint, sep = "_")
  }
  cols <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]
  for (col in cols) {
    if (length(unique(md[[col]])) <= 120) {
      tb <- as.data.frame(table(md[[col]], useNA = "ifany"), stringsAsFactors = FALSE)
      names(tb) <- c("value", "n_cells")
      tb$field <- col
      tb$object <- label
      write.csv(tb, file.path("results", paste0("external_", label, "_", col, "_counts.csv")), row.names = FALSE, fileEncoding = "UTF-8")
      cat("\nFIELD", col, "\n"); print(tb)
    }
  }
  metrics <- names(md)[vapply(md, is.numeric, logical(1)) & grepl("^nFeature|^nCount|percent|mito", names(md), ignore.case = TRUE)]
  grp <- intersect(c("sample", "sample_ID", "orig.ident", "donor", "patient", "timepoint", "donor_timepoint", "disease"), names(md))
  for (col in grp) {
    values <- unique(md[[col]])
    out <- lapply(values, function(g) {
      ix <- if (is.na(g)) is.na(md[[col]]) else !is.na(md[[col]]) & md[[col]] == g
      row <- data.frame(object = label, group_field = col, group_value = as.character(g), n_available_cells = sum(ix))
      for (m in metrics) {
        v <- md[[m]][ix]
        row[[paste0(m, "_median")]] <- median(v, na.rm = TRUE)
        row[[paste0(m, "_min")]] <- min(v, na.rm = TRUE)
        row[[paste0(m, "_max")]] <- max(v, na.rm = TRUE)
      }
      row
    })
    write.csv(do.call(rbind, out), file.path("results", paste0("external_", label, "_qc_available_by_", col, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  }
  rm(obj, md); gc()
}
summary_file <- if ("--full-vaccine" %in% args) "results/external_full_vaccine_summary.csv" else if ("--sle-final" %in% args) "results/external_SLE_GSE135779_final_summary.csv" else "results/external_object_summary.csv"
write.csv(do.call(rbind, summary_rows), summary_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("\nCURRENT AUDIT ENVIRONMENT; NOT HISTORICAL SOFTWARE VERSIONS\n")
print(sessionInfo())
