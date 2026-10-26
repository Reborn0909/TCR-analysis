# TCR-analysis

Reproducible code, data and figures for the single-cell TCR and module-score
analyses: TCR repertoire diversity, donor-level signature scoring, and the
main and supplementary figure sets.

Everything in `analysis/` and `main_figure/` is the original analysis code
moved here unchanged. Comments were removed and absolute data paths replaced
by `DATA_ROOT`; the computational logic is untouched.

## Layout

```
analysis/       26 original R scripts, the analysis pipeline
main_figure/     5 original scripts for the main figure
data/            result tables and signature gene sets
figures/         figures
original/        source inventory (signature gene lists)
environment/     frozen R package inventory
```

## Environment

| Component | Version |
|---|---|
| R | 4.5.2 |
| Seurat | 5.5.0 |
| SeuratObject | 5.3.0 |
| ggplot2 | 4.0.3 |
| ggpubr | 0.6.3 |
| edgeR | Bioconductor |
| Matrix | bundled |

Full inventory: `environment/R_package_inventory.csv`.

## Input data

Single-cell RNA-seq plus matched TCR data from six samples
(H-CD8Trp-1/2, E-CD8Trp-1/2, E-CD8DP-1/2): 76,366 cells x 34,454 genes,
group sizes H-Trp 23,896 / E-Trp 25,613 / E-DP 26,857.

The Seurat objects are 2.5 GB each and are not distributed here. Set
`DATA_ROOT` to the directory holding them before running any script.

## Analysis scripts

| Script | Purpose |
|---|---|
| `01_probe_legacy.R` | inspect legacy object metadata |
| `01_probe_main_inputs.R` | inspect main input objects |
| `02_main_donor_reanalysis.R` | donor-level expression and scoring |
| `02_finalize_from_aggregates.R` | assemble the final result tables |
| `03_current_software_audit.R` | record the software environment |
| `03_external_object_audit.R` | audit external objects |
| `04_depth_matched_TCR.R` | equal-cell TCR resampling, Shannon diversity |
| `09_recover_preQC.R` | recover the pre-QC object |
| `10_recompute_QC_from_raw.R` | recompute QC from raw counts |
| `11_doublet_assessment.R` | doublet assessment |
| `20_panels_changed_new.R` | panel construction |
| `21_panels_english.R` | English panel set |
| `22_volcano_redraw.R` | volcano plots |
| `23_exhaustion_umap_feature.R` | exhaustion UMAP and FeaturePlot |
| `24_exhaustion_marker_panel.R` | exhaustion marker panel |
| `25_exhaustion_marker_figures.R` | exhaustion marker figures |
| `26_featureplot_yellow.R` | FeaturePlot colour scheme |
| `28_heatmap_morandi.R` | Morandi heatmap |
| `29_teacher_marker_panels.R` | marker panels G1 to G6 |
| `30_teacher_panels_extended.R` | extended panels G7 and G8 |
| `32_verify_panel_numbers.R` | verify panel numbers |
| `33_teacher_panels_umap_dotplot.R` | DotPlot, UMAP and violin panels |
| `33b_fix_dotplot.R` | DotPlot fixes |
| `33c_fix_umap_titles.R` | UMAP title fixes |
| `34b_probe_mixing.R` | probe mixing diagnostics |
| `zz_repro_g4.R` | panel G4 reproduction |

## Main figure scripts

`main_figure/01_extract_fast.R`, `01_extract_scores.R`, `02_violin_main_figure.R`
(R) and `02b_violin_local.py`, `03_violin_I_style.py` (Python, matplotlib)
build the six-signature violin figure from the per-cell module scores.

## Signature definitions

`data/signature_gene_sets.csv` — the six main-figure signatures.
`data/signature_gene_sets_ten.csv` — all ten signatures, extracted from
`original/Table_Signature_Gene_Lists.xlsx`.

The CD8 Trp signature is the five-gene core set HLA-DRA, HLA-DRB1, TIGIT,
LAG3, CD74 (Guo et al., Advanced Science 2024).

Scoring: `Seurat::AddModuleScore(assay = "RNA", nbin = 24, ctrl = 100 ,
seed = 20261004)`.

## TCR analysis

`analysis/04_depth_matched_TCR.R` performs equal-cell resampling of the TCR
repertoire: each sample is downsampled to a common cell number, Shannon
diversity is recomputed 2,000 times, and the 2.5th to 97.5th sampling
percentiles are reported. Output is in `data/TCR_equal_cell_downsample_summary.csv`
and `data/TCR_donor_trends.csv`; the figure is
`figures/TCR_equal_cell_downsample.pdf`.

Both donors show the same direction of change at every VDJ depth floor, so
the diversity increase is not explained by differences in recovered TCR depth.
The bars are sampling percentiles from the downsampling procedure, not donor
confidence intervals; with two donors no donor-level interval is estimable.

## Statistics

Two donors. The donor is the biological replicate, so all inference is
donor-level. Cell-level P values printed on the violin plots are descriptive
and are not used for any claim of donor-level significance.

## Integrity

```bash
sha256sum -c MANIFEST_SHA256.csv
```
