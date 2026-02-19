# Shell_scripts_to_generate_features_from_bams

## Purpose
Generate fragmentomic features directly from BAM files in an HPC environment.

## Subfolders
- `1_Runner_scripts/`: user-facing job script generators/launchers.
- `2_Files_scripts_call/`: helper scripts and reference resources called by runner scripts.

## Inputs
- BAM files, reference genome FASTA, peak sets, and filtering resources.

## Outputs
- DELFI ratio outputs, end motif outputs, nucleosome outputs, insert-size summaries.

## Connection
This is the raw feature generation stage for the fragmentomics branch.
