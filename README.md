## cfMeDIP-seq_Data_Resource_Codes

## About
This repository houses the scripts used in the study "A pan-cancer compendium of 1,294 plasma cell-free DNA methylomes and fragmentomes." All data locations referenced by these scripts can be found in Supplementary Table 14 of the manuscript.

<img width="904" height="333" alt="image" src="https://github.com/user-attachments/assets/eb46c93b-d740-4103-a463-5bd86547351e" />

## Directory overview
- **1_Methylation_Scripts/** – generation and processing of methylation-based features
  - `R_scripts_to_process_features/` – numbered R scripts for QC, feature aggregation, batch correction, PCA/UMAP visualization and differential methylation analyses
  - `Shell_scripts_to_generate_features/` – sbatch script to compute methylation features from raw data on a cluster
- **2_Fragmentomics_Scripts/** – extraction of fragmentomic metrics from sequencing BAM files
  - `R_scripts_to_process_features/` – R scripts for fragment length proportions, end motif analysis, nucleosome peak distances and integration of metrics (some scripts have validation cohort variants)
  - `Shell_scripts_to_generate_features_from_bams/`
    - `1_Runner_scripts/` – sbatch runners that launch the pipeline
    - `2_Files_scripts_call/` – helper scripts and references used by runners (includes subfolders `end_motif/`, `nucleosome_peak/`, and `ratio/` with its own README)
- **3_Machine_Learning_Scripts/** – pipelines for building and evaluating classifiers
  - `1_Select_and_PCA_transform_features/` – selection of training/validation sets and PCA transformations
  - `2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/`
    - `1_Runner_scripts/` – shell scripts that initiate full ML pipelines
    - `Scripts_which_runners_call/` – configuration files and R scripts implementing the models
  - `3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/` – Python scripts and CSV configs for cancer type and subtype classification
  - `4_Machine_learning_plotting_scripts/` – R scripts that aggregate results and generate figures

## License
Distributed under the MIT License. See `LICENSE` for more information.

## Contact
Yong Zeng and Dory Abelman

Emails: yzeng@uhn.ca and Dory.Abelman@uhn.ca
