SHELL := /bin/bash
.DEFAULT_GOAL := help

help:
	@echo "Targets:"
	@echo "  make setup      - Install Python requirements and print renv bootstrap commands"
	@echo "  make validate   - Basic repo/documentation sanity checks"
	@echo "  make run        - Print canonical run order and key commands"
	@echo "  make run-cn     - Run cancer-vs-normal script generator (requires config/path setup)"
	@echo "  make run-origin - Run one cancer-origin classification example"
	@echo "  make run-subtype- Run one subtype classification example"

setup:
	python3 -m venv .venv
	. .venv/bin/activate && pip install -r requirements.txt
	@echo "R setup (manual):"
	@echo "  R -q -e 'install.packages(\"renv\"); renv::init(bare = TRUE)'"

validate:
	test -f README.md
	test -f docs/PIPELINE_OVERVIEW.md
	test -f docs/INPUTS.md
	test -f docs/OUTPUTS.md
	test -f docs/REPRODUCIBILITY.md
	test -f config/pipeline_paths.example.env
	test -f requirements.txt
	@echo "Validation checks passed."

run:
	@echo "Ordered run steps:"
	@echo "1) Methylation feature generation shell runner"
	@echo "2) Methylation numbered R scripts"
	@echo "3) Fragmentomics shell runner(s)"
	@echo "4) Fragmentomics numbered R scripts"
	@echo "5) PCA prep scripts (6.2A/6.2B)"
	@echo "6) CN classifier runners"
	@echo "7) Cancer/subtype Python classifiers"
	@echo "8) Plotting scripts"

run-cn:
	bash 3_Machine_Learning_Scripts/2_Shell_and_R_scripts_for_running_ML_pipelines_PE_data_cancer_vs_normal_classifier/1_Runner_scripts/Run_CN_classifier.sh

run-origin:
	@if [[ -z "$$CFMEDIP_MAIN_DIR" ]]; then echo "Set CFMEDIP_MAIN_DIR first."; exit 1; fi
	python 3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/cancer_origin_clf.py motif "Breast Cancer"

run-subtype:
	@if [[ -z "$$CFMEDIP_MAIN_DIR" ]]; then echo "Set CFMEDIP_MAIN_DIR first."; exit 1; fi
	python 3_Machine_Learning_Scripts/3_Scripts_for_running_ML_pipelines_PE_data_cancer_type_or_subtype/subtype_clf.py "IDH mutant Glioma" "IDH mutant Glioma;IDH wildtype Glioma" motif
