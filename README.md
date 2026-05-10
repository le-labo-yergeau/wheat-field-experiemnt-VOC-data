# wheat-field-experiemnt-VOC-data
Reduced Rainfall Impacts on Wheat Rhizosphere Microbiome and Volatile Metabolites

This repository contains the bioinformatics and statistical workflow used to analyze volatile organic compounds (VOCs) and metabolomic responses in the wheat rhizosphere under reduced rainfall conditions. The project focuses on differential abundance analysis, treatment-specific VOC responses, and visualization of drought-associated metabolite patterns across wheat genotypes.
The workflow includes:
Differential abundance analysis using DESeq2
Statistical comparison across rainfall treatments
Identification of significantly produced or consumed VOCs
Generation of publication-quality figures
Creation of master summary tables for downstream interpretation
Experimental Design
Two wheat genotypes were evaluated under different rainfall treatments:
DS — drought-sensitive genotype
DT — drought-tolerant genotype
Rainfall treatments:
100% rainfall (control)
75% rainfall
50% rainfall
25% rainfall
Differential analyses were performed relative to the 100% control condition.
Workflow Summary
1. Differential Abundance Analysis
Differential analysis was performed using the R package:
DESeq2
Contrasts included:
Comparison	Description
DS: 25 vs 100	Severe drought effect in DS
DS: 50 vs 100	Moderate drought effect in DS
DS: 75 vs 100	Mild drought effect in DS
DT: 25 vs 100	Severe drought effect in DT
DT: 50 vs 100	Moderate drought effect in DT
DT: 75 vs 100	Mild drought effect in DT
DT vs DS	Genotype comparison
2. Statistical Outputs
The workflow generates:
Ordered DESeq2 result tables
Significant feature summaries
Master long-format tables
Split metabolite annotations
SVG publication figures
