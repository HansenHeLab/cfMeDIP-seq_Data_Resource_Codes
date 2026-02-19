# 1_Fragmentomic_feature_processing_and_plotting_scripts

## Purpose
Primary cohort processing pipeline for major fragmentomic feature classes.

## Key Scripts
1. `1.1_Fragment_proportions_get_sample_medians.R`
2. `1.2_Fragment_proportions_process_and_plot.R`
3. `2.1_Process_and_plot_insert_sizes.R`
4. `2.2_NMF_analysis_on_insert_sizes.R`
5. `3.1_Process_and_plot_end_motifs.R`
6. `4.1_Process_and_plot_nucelosome_peak_distances.R`
7. `5.1_Integrate_and_plot_zscores_from_different_fragmentomic_features.R`
8. `5.2_Getting_additional_stats_on_integrated_fragmentomic_data.R`

## Inputs
- Fragmentomics matrices (ratios, insert-size, motifs, nucleosome features) and metadata.

## Outputs
- Processed tables, integrated feature scores, and figure panels.

## Connection
Feeds feature summaries used in downstream ML matrix assembly.
