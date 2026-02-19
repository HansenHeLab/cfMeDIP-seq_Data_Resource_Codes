# 3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype

## Purpose
Python workflows for one-vs-all cancer-type and subtype classification.

## Key Files
- `cancer_origin_clf.py`: one-vs-all cancer-type models.
- `subtype_clf.py`: one-vs-all subtype models.
- `run_tab_cancers.csv`: batch combinations for cancer-type runs.
- `run_tab_subtypes.csv`: batch combinations for subtype runs.

## Inputs
- `${CFMEDIP_MAIN_DIR}/data/metadata_df.csv`
- `${CFMEDIP_MAIN_DIR}/data/motif_merged.csv`
- `${CFMEDIP_MAIN_DIR}/data/methyl_merged.csv`
- `${CFMEDIP_MAIN_DIR}/data/delfi_merged.csv`
- `${CFMEDIP_MAIN_DIR}/data/insert_size_merged.csv`
- `${CFMEDIP_MAIN_DIR}/data/ns_peaks_merged.csv`

## Outputs
- `${CFMEDIP_MAIN_DIR}/results_cancers/*.sum.csv` and `*.roc.csv`
- `${CFMEDIP_MAIN_DIR}/results_subtypes/*.sum.csv` and `*.roc.csv`

## Configuration
- `CFMEDIP_MAIN_DIR`: required path root.
- `CFMEDIP_MODELS`: optional model subset for faster tests.
