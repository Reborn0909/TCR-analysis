# What runs, what does not

Every script in this repository was executed against the data that ships here
or against the parent Seurat objects, and the outcome is recorded below.
Scripts are grouped by what they need.

## Environment variables

Two variables control where the scripts look for data:

```bash
export DATA_ROOT=/path/to/dir/holding/seurat_filtered_newQC.rds
```

Scripts read `Sys.getenv("DATA_ROOT", ".")`, so with the variable unset they
resolve paths relative to the repository root.

## Group 1 - runs with the data in this repository

These need only `results/*.csv`, which are included. Verified end to end;
each produced its figure files.

| Script | Output |
|---|---|
| `analysis/21_panels_english.R` | `figures/panels_en/P01` to `P04` |
| `analysis/25_exhaustion_marker_figures.R` | `figures/exhaustion/E1` to `E3` |
| `analysis/22_volcano_redraw.R` | volcano panels |

Run with:

```bash
mkdir -p figures/panels_en figures/exhaustion figures/panels
Rscript analysis/21_panels_english.R
Rscript analysis/25_exhaustion_marker_figures.R
Rscript analysis/22_volcano_redraw.R
```

Output directories must exist; the scripts do not create them.

## Group 2 - needs the parent Seurat object

These call `readRDS(file.path(Sys.getenv("DATA_ROOT", "."), "seurat_filtered_newQC.rds"))`.
That object is 2.48 GB and is not distributed here.

| Script |
|---|
| `analysis/09_recover_preQC.R` |
| `analysis/23_exhaustion_umap_feature.R` |
| `analysis/24_exhaustion_marker_panel.R` |
| `analysis/26_featureplot_yellow.R` |
| `analysis/29_teacher_marker_panels.R` |
| `analysis/30_teacher_panels_extended.R` |
| `analysis/32_verify_panel_numbers.R` |
| `analysis/33_teacher_panels_umap_dotplot.R` |
| `analysis/33b_fix_dotplot.R` |
| `analysis/33c_fix_umap_titles.R` |
| `analysis/34b_probe_mixing.R` |
| `legacy/analysis_script.R` |
| `legacy/NI_Complete_Analysis_Archive.R` |

## Group 3 - needs a file passed on the command line

```bash
Rscript analysis/01_probe_legacy.R <object.rds>
Rscript analysis/01_probe_main_inputs.R <tcr.rds> <expression.rds>
Rscript analysis/02_main_donor_reanalysis.R --input <object.rds> --contig-root <dir> ...
Rscript analysis/03_external_object_audit.R ...
Rscript main_figure/01_extract_fast.R <object.rds> <out.rds>
Rscript main_figure/01_extract_scores.R <object.rds> <out.csv.gz>
```

## Group 4 - platform specific

| Script | Note |
|---|---|
| `analysis/20_panels_changed_new.R` | registers the Windows font `msyh.ttc`; needs that font file or a substitute |
| `analysis/28_heatmap_morandi.R` | needs its expression matrix input |
| `analysis/04_depth_matched_TCR.R` | needs the per-cell contig key tables `cellkeys_*.rds` |
| `legacy/CD8_analysis.R` | reads objects that must already be in memory |
| `legacy/TCR_analysis.R` | reads `seurat_filtered` from the R session; it is a continuation script |

## Not included

Two files from the original collection, `TCR_newQC.R` and `TCR_analysis_plots.R`,
are R console transcripts rather than source files: they contain the `>` and
`+` prompts and interleaved prose, so they do not parse as R code and are not
part of this repository.
