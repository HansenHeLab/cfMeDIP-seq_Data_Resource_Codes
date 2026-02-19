# 2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier

## Purpose
Cancer-vs-normal classifier workflow using shell script generation plus R model execution.

## Subfolders
- `1_Runner_scripts/`: generate job scripts from config CSVs.
- `Scripts_which_runners_call/`: configs and R scripts implementing model training/validation.

## Inputs
- Config CSV files listing feature matrix paths, metadata paths, and output folders.

## Outputs
- Generated job scripts, slurm logs, model performance summaries, and model/prediction objects.

## Connection
This is the binary classification branch feeding downstream result plotting.
