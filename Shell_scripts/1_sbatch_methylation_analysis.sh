#!/bin/bash
#SBATCH -p veryhimem
#SBATCH --mem=800G
#SBATCH -t 5-00:00:00
#SBATCH -J aggr_cnt
#SBATCH -o aggr_cnt_%j.out
#SBATCH -e aggr_cnt_%j.err

##
module load R/4.0.0

cd  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/raw_cnt

## sample aggregation 
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_Normal_SE.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_AML.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_Normal.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_Normal_PE.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_SE.csv
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/aggr_quant.R ../TCGE_cfDM_samples_only_fixed_PE.csv

## batch correction and normalization for healthy samples
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/norm_batch_correct.R

## batch correction and normalization for all SE and PE samples, separately
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/pan_batch_correct_SE_PE.R

## DMR analyses
cd /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs
mkdir -p PE
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/pan_DMR_PE.R

mkdir -p SE
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/pan_DMR_SE.R

mkdir -p PE_Specific
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/cancer_specific_DMR_PE.R

mkdir -p SE_Specific
Rscript --vanilla  /cluster/projects/tcge/cell_free_epigenomics/processed_data/merge/DMRs/cancer_specific_DMR_SE.R
