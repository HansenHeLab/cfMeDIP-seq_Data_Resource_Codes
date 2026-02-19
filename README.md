# cfMeDIP-seq Data Resource Codes

## Overview
This repository contains analysis code used in the study "A pan-cancer compendium of 1,294 plasma cell-free DNA methylomes and fragmentomes." It supports three linked workflows: methylation feature generation/processing, fragmentomics feature generation/processing, and downstream machine-learning analyses (cancer-vs-normal, cancer-type, and subtype classification). The codebase is optimized for HPC-style batch execution and assumes curated input tables from the study data resource.

## Repository Layout
- `1_Methylation_Scripts/`: methylation feature generation (cluster runners) and downstream R analysis scripts.
- `2_Fragmentomics_Scripts/`: fragmentomics feature generation from BAM files and downstream R analysis scripts.
- `3_Machine_Learning_Scripts/`: feature selection/PCA, classifier training, and result plotting.
- `docs/`: end-to-end usage documentation, input/output schema notes, and reproducibility guidance.
- `config/`: local configuration templates for path overrides.
- `Makefile`: lightweight wrappers for setup, validation, and common run commands.

## Quickstart
The repository does not include full study input data. Quickstart therefore means "configure environment + run a minimal command against your own prepared inputs".

```bash
# 1) Clone and enter repo
git clone <repo-url>
cd cfMeDIP-seq_Data_Resource_Codes

# 2) Python environment (for cancer-type/subtype classifiers)
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 3) Optional R environment bootstrap
R -q -e 'install.packages("renv"); renv::init(bare = TRUE)'

# 4) Validate expected files/entrypoints
make validate

# 5) Run one cancer-type classifier job (requires prepared data directory)
export CFMEDIP_MAIN_DIR=/absolute/path/to/ML_workspace
python 3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/cancer_origin_clf.py motif "Breast Cancer"
```

## Installation And Environment Setup

### R (primary analysis language)
Recommended: R 4.1+ on Linux, with explicit package snapshots via `renv`.

```bash
R -q -e 'install.packages("renv")'
R -q -e 'renv::init(bare = TRUE)'
R -q -e 'renv::install(c("tidyverse","dplyr","ggplot2","caret","pROC","doParallel","DESeq2","sva"))'
R -q -e 'renv::snapshot()'
```

Note: many scripts rely on additional Bioconductor and visualization packages. See `docs/REPRODUCIBILITY.md` for expanded package guidance.

### Python (ML subtype/cancer-type scripts)
Use `requirements.txt`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Configuration
Use `config/pipeline_paths.example.env` as a local template:

```bash
cp config/pipeline_paths.example.env config/pipeline_paths.env
# edit values for your system
source config/pipeline_paths.env
```

Supported overrides:
- `CFMEDIP_MAIN_DIR`: root containing `data/` and `results_*` for Python classifiers.
- `CFMEDIP_MODELS`: optional comma-separated model subset for Python scripts (example: `lr,rf`).
- `CFMEDIP_CN_CONFIG`, `CFMEDIP_CN_OUTPUT_DIR`, `CFMEDIP_CN_SLURM_DIR`, `CFMEDIP_CN_RUNNER_SCRIPT`, `CFMEDIP_CN_SCRIPTS_DIR`: overrides for cancer-vs-normal shell/R runners.

## Running The Full Pipeline (Ordered)
1. Generate methylation features from raw inputs:
   - `1_Methylation_Scripts/Shell_scripts_to_generate_features/1_sbatch_methylation_analysis.sh`
2. Process methylation features in numbered order (`1_...` through `5.6_...`):
   - `1_Methylation_Scripts/R_scripts_to_process_features/`
3. Generate fragmentomics features from BAMs (choose runner for your cohort):
   - `2_Fragmentomics_Scripts/Shell_scripts_to_generate_features_from_bams/1_Runner_scripts/`
4. Process fragmentomics features in numbered order (`1.x` to `5.x`, then optional `9.x/10.x` statistics):
   - `2_Fragmentomics_Scripts/R_scripts_to_process_features/`
5. Build train/validation feature matrices and PCA transforms:
   - `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/6.2A_...PE...R`
   - `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/6.2B_...SE...R`
6. Run cancer-vs-normal classifiers:
   - `3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/1_Runner_scripts/Run_CN_classifier.sh`
   - `.../Run_CN_classifier_SE.sh`
7. Run cancer-type and subtype classifiers (Python):
   - `3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/cancer_origin_clf.py`
   - `3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/subtype_clf.py`
8. Generate downstream summary plots:
   - `3_Machine_Learning_Scripts/4_Machine_learning_plotting_scripts/`

Detailed workflow text and flowchart are in `docs/PIPELINE_OVERVIEW.md`.

## Outputs
Primary outputs include:
- Processed feature matrices (`.rds`, `.csv`) from methylation/fragmentomics processing.
- PCA-transformed training/validation matrices for ML.
- Classifier outputs (`*.sum.csv`, `*.roc.csv`, model objects `.rds`, kappa summaries `.txt/.pdf`).
- Publication-style figure panels from plotting scripts.

See `docs/OUTPUTS.md` for file naming conventions and interpretation.

## Troubleshooting
- `File not found` errors: most scripts contain hardcoded absolute paths from the original compute environment; replace these with your local paths before running.
- Empty intersections of metadata/features: ensure sample IDs are normalized consistently (for example remove `_dedup` where expected).
- Cross-validation failures (`n_splits` too large): class sizes are too small for requested folds. Reduce cohort scope or use `CFMEDIP_MODELS` to run simpler checks first.
- R package missing: install package into your `renv` or system library and rerun.
- HPC module errors: cluster runner scripts assume `module load` and scheduler directives are available.

## Additional Documentation
- `docs/PIPELINE_OVERVIEW.md`
- `docs/INPUTS.md`
- `docs/OUTPUTS.md`
- `docs/REPRODUCIBILITY.md`
- `docs/DOCS_IMPROVEMENTS_SUMMARY.md`
