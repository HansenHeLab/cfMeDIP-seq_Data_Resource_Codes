# R_scripts_to_process_features

## Purpose
Post-processing and statistical analysis of fragmentomic features generated from BAM-derived outputs.

## Subfolders
- `1_Fragmentomic_feature_processing_and_plotting_scripts/`: primary cohort processing.
- `2_Fragmentomic_feature_processing_and_plotting_scripts_validation_cohort/`: validation cohort processing.
- `3_Other_scripts_for_calculating_statistics/`: additional revision/statistical analyses.

## Inputs
- Fragmentomic matrices from shell feature generation outputs.
- Metadata and cohort labels referenced in scripts.

## Outputs
- Feature-level summaries, integrated z-score matrices, plots, and cohort statistics.

## Upstream/Downstream
- Upstream: `../Shell_scripts_to_generate_features_from_bams/1_Runner_scripts/`.
- Downstream: ML preparation scripts in `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/`.
