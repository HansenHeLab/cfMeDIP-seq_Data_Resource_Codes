# 3_Primary_runner_scripts_scriptgens_call

## Purpose
Primary orchestrator script for CN model benchmarking/training.

## Key Script
- `Runner_fragmentation_Cancer_vs_normal.R`

## Inputs
- Command-line arguments from generated shell scripts (technology name, data paths, labels, sample ID column, optional external validation path).

## Outputs
- Kappa summary plot/table, CV predictions, fitted model objects, optional external validation outputs.

## Configuration
Supports `CFMEDIP_CN_SCRIPTS_DIR` to locate secondary helper scripts.
