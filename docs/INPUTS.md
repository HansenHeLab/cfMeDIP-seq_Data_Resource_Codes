# Inputs

## Input Classes
1. Raw sequencing data (BAM and count-level inputs for feature generation).
2. Preprocessed feature matrices (`.rds` in R pipelines, `.csv` in Python pipelines).
3. Metadata tables with consistent sample identifiers.
4. Run configuration tables (`.csv`) for batch model execution.

## Required Metadata Fields

### Cancer-Type Python Classifier (`cancer_origin_clf.py`)
`metadata_df.csv` must include:
- `sample_id`
- `cancer_type`
- `project_id`
- `type`

Behavior:
- Excludes rows where `cancer_type == "Normal"`.
- Excludes rows where `project_id == "TCGE-CFMe-HCC"`.
- Keeps only `type == "PE"`.
- Removes `_dedup` suffix from `sample_id` before matching.

### Subtype Python Classifier (`subtype_clf.py`)
`metadata_df.csv` must include:
- `sample_id`
- `cancer_subtype`
- `type`

Behavior:
- Keeps only subtype values supplied via CLI argument.
- Keeps only `type == "PE"`.
- Removes `_dedup` suffix from `sample_id` before matching.

### Cancer-vs-Normal R Runner (`Runner_fragmentation_Cancer_vs_normal.R`)
Metadata RDS is expected to include:
- `CN_classifier` (binary group labels such as `healthy`/`cancer`)
- a sample ID column matching `sample_id_column` passed from config (for example `sample_id`, `endmotif_name`, `sequencing_id`)

## Feature Matrix Requirements
- Rows are samples; columns are features.
- Row names (RDS) or index column (CSV) must be sample IDs.
- No duplicate sample IDs.
- For merged feature runs (`all`, `motif+methyl`, etc.), scripts use inner joins; only common samples survive.

## Run Tables And Config Files

### `run_tab_cancers.csv`
Required columns:
- `Feature`
- `Cancer`

### `run_tab_subtypes.csv`
Required columns:
- `features`
- `target_subtype`
- `all_subtypes_in_group` (semicolon-separated)
- `group`

### `Machine_learning_config_combined_CN_classifier_*.csv`
Required columns:
- `technology_name`
- `input_data`
- `metadata`
- `output_folder`
- `title`
- `title2`
- `sample_id_column`
- `external_data_file`

## Naming And ID Conventions
- Keep sample IDs stable across all modalities.
- Apply any suffix normalization (such as `_dedup`) consistently before merging.
- Ensure label values match script expectations exactly (`healthy`, `cancer`, cancer/subtype strings).
