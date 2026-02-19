# Shell_scripts_to_generate_features

## Purpose
Cluster submission script for methylation feature aggregation and first-pass preprocessing.

## Key Script
- `1_sbatch_methylation_analysis.sh`: orchestrates aggregation, normalization/batch correction, and DMR script entrypoints in the original HPC environment.

## Inputs
- Count files and sample list CSVs referenced in script path variables.

## Outputs
- Aggregated count matrices and DMR-ready intermediate files.

## Workflow Position
Run first in methylation workflow, before numbered R analysis scripts.
