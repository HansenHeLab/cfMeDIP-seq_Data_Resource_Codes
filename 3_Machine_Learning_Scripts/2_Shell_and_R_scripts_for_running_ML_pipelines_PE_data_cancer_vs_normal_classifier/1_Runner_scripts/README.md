# 1_Runner_scripts

## Purpose
Generate cancer-vs-normal cluster scripts from config tables.

## Key Scripts
- `Run_CN_classifier.sh`
- `Run_CN_classifier_SE.sh`

## Inputs
- `Machine_learning_config_combined_CN_classifier_*.csv` tables.

## Outputs
- One generated shell script per technology/group comparison.

## Configuration
Supports path overrides via environment variables:
- `CFMEDIP_CN_CONFIG`, `CFMEDIP_CN_OUTPUT_DIR`, `CFMEDIP_CN_SLURM_DIR`, `CFMEDIP_CN_RUNNER_SCRIPT`, `CFMEDIP_R_MODULE`.
