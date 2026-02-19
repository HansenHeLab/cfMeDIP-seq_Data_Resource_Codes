# R_scripts_to_process_features

## Purpose
Numbered methylation processing workflow from sample metadata/QC through DMR calls and filtering.

## Suggested Run Order
1. `1_sample_info_summary.R`
2. `2_QC.R`
3. `3_Aggregate_quantifications.R`
4. `4.1_Batch_correction_for_healthy.R`
5. `4.2_Batch_correction_for_SE_PE.R`
6. `4.3_PCA_UMAP_plots.R`
7. `5.1_pan-cancer_DMR.R` and `5.1.1_pan-cancer_DMR_age_sex.R`
8. `5.2_cancer_specific_DMR.R`
9. `5.3_DMR_filtering_blood_age_sex.R` and `5.3.1_DMR_filtering_blood_age_sex_assess.R`
10. `5.4.1_annotations_and_permutaion_functions.R` and `5.4_DMR_enrichment.R`
11. `5.5_pan_cancer_DMR_after_filtering.R`
12. `5.6_cancer_sepcific_DMR_after_filtering.R`

## Inputs
- Intermediate methylation count matrices (`.rds`, `.csv`) and metadata files referenced in each script.

## Outputs
- QC plots/tables, normalized matrices, DMR result objects, and enrichment outputs.

## Upstream/Downstream
- Upstream: `../Shell_scripts_to_generate_features/1_sbatch_methylation_analysis.sh`.
- Downstream: PCA/ML scripts in `3_Machine_Learning_Scripts/1_Select_and_PCA_transform_features/`.
