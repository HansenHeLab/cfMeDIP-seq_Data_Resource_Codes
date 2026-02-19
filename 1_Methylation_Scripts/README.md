# 1_Methylation_Scripts

## Purpose
Scripts for generating methylation feature counts and performing downstream methylation analyses (QC, batch correction, DMR discovery, and enrichment).

## Subfolders
- `R_scripts_to_process_features/`: numbered R workflow scripts for methylation feature processing and DMR analyses.
- `Shell_scripts_to_generate_features/`: cluster runner for initial count aggregation and core preprocessing.

## Upstream Inputs
- Raw count inputs and sample list CSV files from the cfMeDIP-seq data resource.
- Metadata tables used by downstream R scripts.

## Downstream Outputs
- Batch-corrected and normalized methylation matrices.
- DMR results and enrichment summaries used by ML feature-prep scripts.

## Workflow Connection
Run this module before `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/` when methylation features are needed for classification.
