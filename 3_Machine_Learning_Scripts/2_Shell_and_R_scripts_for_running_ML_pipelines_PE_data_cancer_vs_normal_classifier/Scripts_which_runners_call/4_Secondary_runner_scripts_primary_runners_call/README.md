# 4_Secondary_runner_scripts_primary_runners_call

## Purpose
Secondary helper scripts sourced by the primary CN runner.

## Files
- `FuncClassifier.R`: utility splitting/downsampling helpers.
- `classifier_opt.R`: model optimization helper logic.
- `final_classifier.R`: final classifier fit/CV/external-validation workflow.

## Connection
These scripts are sourced by `Runner_fragmentation_Cancer_vs_normal.R` and are not typically invoked standalone.
