# Outputs

## Overview
Outputs are generated across workflow stages and are mostly written as `.rds`, `.csv`, `.txt`, and `.pdf` files. Output roots are often controlled by hardcoded script variables or path overrides.

## Key Output Groups

### 1) Methylation Processing Outputs
Typical artifacts:
- QC summaries and diagnostic plots.
- Batch-corrected / normalized matrices.
- DMR result objects and enrichment summaries.

Primary locations:
- Directories defined inside numbered R scripts under `1_Methylation_Scripts/R_scripts_to_process_features/`.

### 2) Fragmentomics Processing Outputs
Typical artifacts:
- Fragment ratio matrices.
- Insert-size/end-motif/nucleosome summaries.
- Integrated z-score matrices and statistics plots.

Primary locations:
- Directories defined inside R scripts under `2_Fragmentomics_Scripts/R_scripts_to_process_features/`.
- For raw feature extraction: output roots declared in shell runners under `2_Fragmentomics_Scripts/Shell_scripts_to_generate_features_from_bams/1_Runner_scripts/`.

### 3) PCA/Training Matrix Outputs
From `6.2A` and `6.2B` scripts:
- PCA objects (`*_PCA_Object.rds`)
- Training/validation PCA score matrices (`*_Training_PCA_*percent.rds`, `*_Validation_PCA_*percent.rds`)
- Combined multi-feature PCA matrices.
- Scree and cumulative variance plots.

### 4) Cancer-Vs-Normal Classifier Outputs
From `Runner_fragmentation_Cancer_vs_normal.R`:
- Model comparison table: `*_results_<type1>_vs_<type2>.txt`
- Kappa plot: `*_kappa_<type1>_vs_<type2>.pdf`
- Internal CV predictions: `*_outcomes_<type1>_vs_<type2>.rds`
- Optional external validation predictions and final model objects.

### 5) Cancer-Type/Subtype Python Outputs
From `cancer_origin_clf.py` and `subtype_clf.py`:
- Summary table: `<label>+<feature>.sum.csv`
- ROC trace table: `<label>+<feature>.roc.csv`

Default target directories:
- `${CFMEDIP_MAIN_DIR}/results_cancers/`
- `${CFMEDIP_MAIN_DIR}/results_subtypes/`

## Output Interpretation Notes
- `.sum.csv` files are fold-level model performance summaries.
- `.roc.csv` files provide points for ROC plotting (FPR/TPR per run).
- `.rds` model objects can be reloaded for downstream diagnostics and figure regeneration.
