# 2_Fragmentomics_Scripts

## Purpose
Feature extraction and analysis code for cfDNA fragmentomics (ratios, insert sizes, end motifs, nucleosome distances, integrated signatures).

## Subfolders
- `Shell_scripts_to_generate_features_from_bams/`: BAM-level feature generation runners and helper code.
- `R_scripts_to_process_features/`: downstream statistical processing and plotting.

## Upstream Inputs
- BAM files and reference resources (genome FASTA, peak sets, filtered regions).

## Downstream Outputs
- Processed fragmentomic matrices and integrated feature summaries for machine learning.

## Workflow Connection
Run this module before ML feature preparation in `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/`.
