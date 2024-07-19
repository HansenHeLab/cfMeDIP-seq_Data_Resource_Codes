#!/bin/bash

# Define cancer types
declare -a types=("normal" "metastatic_castration_resistant_prostate_cancer" "head_and_neck_squamous_cell_carcinoma" "meningioma" "low-grade_glioneuronal" "idh_wildtype_glioma" "idh_mutant_glioma" "brain_metastases" "hemangiopericytoma" "aml" "primary_prostate_cancer" "UM" "UM_Relapse" "lfs_survivor" "lfs_previvor" "lfs_positive" "samll_cell_lung_cancer")

# Paths
config_file="/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/June_2024/Machine_learning_config_end_methylation_cancer_subtype.csv"
output_dir="/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/June_2024/output_scripts_cancer_subtype_combined"
slurm_dir="/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/June_2024/slurm_outputs_cancer_subtype_combined"

# Create directories if they do not exist
mkdir -p "$output_dir"
mkdir -p "$slurm_dir"

# Function to trim whitespace and special characters
trim() {
    local var="$*"
    # Remove leading whitespace
    var="${var#"${var%%[![:space:]]*}"}"
    # Remove trailing whitespace
    var="${var%"${var##*[![:space:]]}"}"
    echo -n "$var"
}

# Read the config file and generate SLURM scripts
{
  read
  while IFS=, read -r technology_name input_data metadata output_folder title title2 sample_id_column
  do
    technology_name=$(trim "$technology_name")
    input_data=$(trim "$input_data")
    metadata=$(trim "$metadata")
    output_folder=$(trim "$output_folder")
    title=$(trim "$title")
    title2=$(trim "$title2")
    sample_id_column=$(trim "$sample_id_column")

    for i in "${!types[@]}"; do
      for j in $(seq $((i + 1)) $((${#types[@]} - 1))); do
        type1="${types[i]}"
        type2="${types[j]}"

        script_file="$output_dir/${technology_name}_${type1}_vs_${type2}.sh"
        
        cat <<EOL > "$script_file"
#!/bin/bash
#SBATCH --job-name=${technology_name}_${type1}_vs_${type2}
#SBATCH --output=${slurm_dir}/${technology_name}_${type1}_vs_${type2}.out
#SBATCH --error=${slurm_dir}/${technology_name}_${type1}_vs_${type2}.err
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --mem=8G

module load R/4.1

Rscript /cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run2/Machine_Learning/June_2024/Runner_fragmentation_cancer_subtype.R "${technology_name}" "${input_data}" "${metadata}" "${output_folder}" "${title}" "${title2}" "${sample_id_column}" "${type1}" "${type2}"
EOL
      done
    done
  done
} < "$config_file"

