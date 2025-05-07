#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

INPUTDIR=/cluster/projects/tcge/cell_free_epigenomics/processed_data/TCGE-CFMe-MCA/dedup_bam_se
base=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run5_single_end_data
outdir_endmotif=$base/Output_Endmotif
ref=/cluster/projects/tcge/DB/MEDIPIPE/hg38/bwa_index_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
shdir=$base/sh_scripts_endmotif
slurm=$base/slurm_endmotif
updated_motif_script_dir=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory
fragmentomics_endmotif=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/updated_endmotif_scripts

mkdir -p $base
mkdir -p $outdir_endmotif
mkdir -p $shdir
mkdir -p $slurm

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}_endmotif.sh
#!/bin/bash
#
#$ -cwd

#SBATCH --job-name=${bam}_endmotif
#SBATCH --output=$slurm/${bam}_endmotif.out
#SBATCH --error=$slurm/${bam}_endmotif.err
#SBATCH --time=2-12:00:00
#SBATCH --mem=20GB
#SBATCH -c 1

module load samtools
module load bedtools/2.27.1
module load R/4.0.0

#### Convert BAM to BED format for all reads
samtools view -b $INPUTDIR/${bam}.bam | \
bedtools bamtobed -i stdin > $outdir_endmotif/${bam}.bed

## Check if BED file is not empty
if [ -s $outdir_endmotif/${bam}.bed ]; then

    ## Format BED for end motif analysis
    Rscript $updated_motif_script_dir/motif_get_bed.R \
    --id ${bam} \
    --bed $outdir_endmotif/${bam}.bed \
    --outdir $outdir_endmotif

    ## Get FASTA sequences
    bedtools getfasta -bedOut -fi $ref -bed $outdir_endmotif/${bam}_5.bed > $outdir_endmotif/${bam}_fasta_5.bed
    bedtools getfasta -bedOut -fi $ref -bed $outdir_endmotif/${bam}_3.bed > $outdir_endmotif/${bam}_fasta_3.bed

    ## Convert FASTA to end motif context frequencies
    Rscript $fragmentomics_endmotif/motif_get_contexts.R \
    --id ${bam} \
    --fasta_5 $outdir_endmotif/${bam}_fasta_5.bed \
    --fasta_3 $outdir_endmotif/${bam}_fasta_3.bed \
    --outdir $outdir_endmotif

    ## Clean up intermediate files
    if [[ -f "$outdir_endmotif/${bam}_raw.txt" ]]
    then
        echo -e "Endmotif output completed successfully for $bam"
        rm $outdir_endmotif/${bam}*.bed*
    else
        echo -e "Errors in end motif run for $bam"
    fi

else
    echo "No reads found for $bam"
fi

EOF
done
