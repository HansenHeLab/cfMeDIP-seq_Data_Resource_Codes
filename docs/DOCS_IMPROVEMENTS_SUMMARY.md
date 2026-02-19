# Documentation Improvements Summary

## Scope Completed
This pass focused on usability and documentation, with additive changes only. Scientific analysis logic was not changed. Small usability-only path/config overrides were added to selected runner scripts.

## Files Changed

### Top-level and docs
- `README.md`
- `Makefile`
- `requirements.txt`
- `config/pipeline_paths.example.env`
- `docs/PIPELINE_OVERVIEW.md`
- `docs/INPUTS.md`
- `docs/OUTPUTS.md`
- `docs/REPRODUCIBILITY.md`
- `docs/DOCS_IMPROVEMENTS_SUMMARY.md`

### Folder-level READMEs added
- `1_Methylation_Scripts/README.md`
- `1_Methylation_Scripts/R_scripts_to_process_features/README.md`
- `1_Methylation_Scripts/Shell_scripts_to_generate_features/README.md`
- `2_Fragmentomics_Scripts/README.md`
- `2_Fragmentomics_Scripts/R_scripts_to_process_features/README.md`
- `2_Fragmentomics_Scripts/R_scripts_to_process_features/1_Fragmentomic_feature_processing_and_plotting_scripts/README.md`
- `2_Fragmentomics_Scripts/R_scripts_to_process_features/2_Fragmentomic_feature_processing_and_plotting_scripts_validation_cohort/README.md`
- `2_Fragmentomics_Scripts/R_scripts_to_process_features/3_Other_scripts_for_calculating_statistics/README.md`
- `2_Fragmentomics_Scripts/Shell_scripts_to_generate_features_from_bams/README.md`
- `2_Fragmentomics_Scripts/Shell_scripts_to_generate_features_from_bams/1_Runner_scripts/README.md`
- `3_Machine_Learning_Scripts/README.md`
- `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/1_Runner_scripts/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/Scripts_which_runners_call/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/Scripts_which_runners_call/2_Configs_scriptgens_call/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/Scripts_which_runners_call/3_Primary_runner_scripts_scriptgens_call/README.md`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/Scripts_which_runners_call/4_Secondary_runner_scripts_primary_runners_call/README.md`
- `3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/README.md`
- `3_Machine_Learning_Scripts/4_Machine_learning_plotting_scripts/README.md`

### Script header documentation updates
Added standardized header comment blocks (Purpose, Inputs, Outputs, How to run, Key parameters, Dependencies) to runnable scripts that were missing them.

Updated script groups include:
- Methylation processing scripts in `1_Methylation_Scripts/R_scripts_to_process_features/`
- Methylation shell runner in `1_Methylation_Scripts/Shell_scripts_to_generate_features/`
- Fragmentomics processing scripts in `2_Fragmentomics_Scripts/R_scripts_to_process_features/**`
- Fragmentomics shell runners in `2_Fragmentomics_Scripts/Shell_scripts_to_generate_features_from_bams/1_Runner_scripts/`
- CN runner/secondary scripts in `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/**`
- Python subtype classifier script and selected ML plotting scripts

### Usability/path configurability updates
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/1_Runner_scripts/Run_CN_classifier.sh`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/1_Runner_scripts/Run_CN_classifier_SE.sh`
- `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/Scripts_which_runners_call/3_Primary_runner_scripts_scriptgens_call/Runner_fragmentation_Cancer_vs_normal.R`
- `3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/cancer_origin_clf.py`
- `3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/subtype_clf.py`

## What Was Added
1. Clear high-level onboarding docs (purpose, structure, run order, setup, troubleshooting).
2. Dedicated docs pages for flow, input schemas, outputs, and reproducibility.
3. Folder-level READMEs explaining local purpose and upstream/downstream links.
4. Per-script header documentation for runnable scripts lacking a header.
5. Lightweight run/validation wrappers via `Makefile`.
6. Explicit configuration template for path overrides.
7. Python dependency file (`requirements.txt`) and path/model override support for ML scripts.

## Ambiguous Areas
1. Many scripts still contain original hardcoded cluster paths and depend on local infrastructure/data not distributed in this repository. Documentation now flags this, but per-site path adaptation remains required.
2. Full manuscript data inputs are not included, so true end-to-end reproduction cannot be executed directly from a clean clone without external data provisioning from the manuscript directory.
