# 1_Select_and_PCA_transform_features

## Purpose
Prepare ML-ready feature matrices by selecting cohorts, splitting train/validation sets, and running PCA transformations.

## Key Scripts
- `6.2A_Selecting_samples_for_training_and_validation_cohorts_PE_data_updatedJan_PCA_transform.R`
- `6.2B_Selecting_samples_for_training_and_validation_cohorts_SE_data_updatedMar_PCA_transform.R`

## Inputs
- Raw/processed feature RDS matrices and metadata RDS files referenced in script variables.

## Outputs
- PCA objects, transformed train/validation matrices, and summary plots/tables.

## Connection
Outputs are consumed by CN, cancer-origin, and subtype classification scripts.
