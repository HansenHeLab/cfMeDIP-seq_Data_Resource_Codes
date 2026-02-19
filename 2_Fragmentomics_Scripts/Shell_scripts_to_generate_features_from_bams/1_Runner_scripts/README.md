# 1_Runner_scripts

## Purpose
HPC entrypoint scripts for generating fragmentomics features from BAMs.

## Key Scripts
- `2_run_sbatch_insert_size_fragment_ratios_and_nucleosome_peak.sh`
- `3_run_sbatch_end_motifs.sh`
- `Run5_Fall2024_cases_HCC.sh`
- `Run5_Fall2024_cases_Inspire.sh`
- `Run5_Fall2024_cases_Inspire_endmotif_updated.sh`
- `Run5_Single_End_Data_endmotif_only.sh`

## Inputs
- BAM directory path configured in each script (`INPUTDIR`).
- Cluster modules/tools (`samtools`, `bedtools`, `picard`, `R`).

## Outputs
- Per-sample fragmentomics output files in runner-defined output directories.

## Connection
Outputs feed into `2_Fragmentomics_Scripts/R_scripts_to_process_features/`.
