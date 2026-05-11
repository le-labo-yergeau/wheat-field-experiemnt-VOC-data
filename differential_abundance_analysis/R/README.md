
# Wheat Rhizosphere VOC Analysis

This repository contains R scripts, processed datasets, figures, and analysis outputs for wheat rhizosphere volatile organic compound (VOC) analysis.

## Files Description

- `VOC.Rmd` — Main R Markdown workflow for VOC statistical analysis and figure generation.
- `VOC.html` — Rendered HTML report generated from the R Markdown analysis.
- `Wheat Field Experiment VOC analysis_f.R` — Standalone R script for VOC data processing and visualization.
- `VOC_analysis.Rproj` — RStudio project file.

## Input Data

- `Field_exp_wheat_alignment0.5_sorted_nist_metabo.xlsx` — Processed metabolomics/VOC annotation dataset.
- `Figure7_ABC_Data.xlsx` — Input data used for Figure 7 visualizations.
- `meta.xlsx` — Sample metadata and experimental design information.

## Output Data

- `heatmap_mean_matrix.csv` — Heatmap abundance matrix used for clustering and visualization.
- `heatmap_mean_matrix_with_metadata.csv` — Heatmap matrix merged with feature metadata and statistics.

## Folders

- `output_figures/` — Generated figures and plots from the VOC analysis workflow.
- `output_tables/` — Statistical result tables and processed outputs.
