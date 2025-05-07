#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

INPUTDIR=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/bam_links_all
base=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run4_allfrags_v2
outdir_delfi=$base/Output_Delfi
shdir=$base/sh_scripts
fragmentomics_delfi=/cluster/projects/pughlab/bin/fragmentomics/v2/ratio # The DELFI plots

ref=/cluster/projects/tcge/DB/MEDIPIPE/hg38/bwa_index_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
picard_dir=/cluster/tools/software/picard/2.10.9
fragmentomics_peak=/cluster/projects/pughlab/bin/fragmentomics/v2/nucleosome_peak
peaks=/cluster/projects/pughlab/bin/fragmentomics/v2/nucleosome_peak/peaks/hg38
outdir_peak=$base/Output_Peak

fragmentomics_endmotif=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/updated_endmotif_scripts
outdir_endmotif=$base/Output_Endmotif

insertsize=$base/Insert_size

slurm=$base/slurm

mkdir -p $base
mkdir -p $outdir_delfi
mkdir -p $outdir_peak
mkdir -p $outdir_endmotif
mkdir -p $shdir
mkdir -p $slurm
mkdir -p $insertsize

module load picard/2.10.9

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}.sh
#!/bin/bash
#
#$ -cwd

#SBATCH --job-name=${bam}_frag
#SBATCH --output=$slurm/${bam}_frag.out
#SBATCH --error=$slurm/${bam}_frag.err
#SBATCH --time=24:00:00
#SBATCH --mem=30GB
#SBATCH -c 1

module load picard/2.10.9
module load igenome-human/hg38
module load samtools
module load R/4.0.0
module load samtools/1.10
module load picard/2.10.9
module load bedtools/2.27.1

## Subset the bam file for insert sizes between 20bp and 600bp
samtools view -h $INPUTDIR/${bam}.bam | awk 'substr(\$0,1,1)=="@" || (\$9>=20 && \$9<=600) || (\$9<=-20 && \$9>=-600)' | \
samtools view -b > $outdir_peak/${bam}.bam
samtools index $outdir_peak/${bam}.bam

## Remove duplicates
java -Xmx15g -jar $picard_dir/picard.jar MarkDuplicates \
I=$outdir_peak/${bam}.bam \
O=$outdir_peak/${bam}_deduped.bam \
M=$outdir_peak/${bam}_metrics.txt \
TMP_DIR=$outdir_peak \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true \
MAX_RECORDS_IN_RAM=100000 

rm $outdir_peak/${bam}.bam

samtools sort -n $outdir_peak/${bam}_deduped.bam -o $outdir_peak/${bam}_deduped_sorted.bam

rm $outdir_peak/${bam}_deduped.bam

## Sort bam by name and convert to bedpe format
samtools view -bf 0x2 $outdir_peak/${bam}_deduped_sorted.bam | \
bedtools bamtobed \
-i stdin \
-bedpe > $outdir_peak/${bam}.bedpe

rm $outdir_peak/${bam}_deduped_sorted.bam

## Split bedpe into chromosomes
#$fragmentomics_peak/R/splitBedpe.sh $outdir_peak/${bam}.bedpe

#### Run end motif - continuing on last steps, 167bp fragments
## Format bedpe to bed
Rscript $fragmentomics_endmotif/motif_format_bedpe.R \
--id $bam \
--bedpe $outdir_peak/${bam}.bedpe \
--outdir $outdir_endmotif

rm $outdir_peak/${bam}.bedpe

## Get FASTA sequences 
bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir_endmotif/${bam}_5.bed > $outdir_endmotif/${bam}_fasta_5.bed
bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir_endmotif/${bam}_3.bed > $outdir_endmotif/${bam}_fasta_3.bed

## Convert FASTA to end motif context frequencies
Rscript $fragmentomics_endmotif/motif_get_contexts.R \
--id $bam \
--fasta_5 $outdir_endmotif/${bam}_fasta_5.bed \
--fasta_3 $outdir_endmotif/${bam}_fasta_3.bed \
--outdir $outdir_endmotif


#### Cleaning up
## Remove intermediate files - end motif
if [[ -f "$outdir_endmotif/${bam}_5prime_raw.txt" ]]
then
echo -e "Endmotif output completed sucessfully"
rm $outdir_peak/${bam}*.bam*
rm $outdir_peak/${bam}_metrics.txt
rm $outdir_peak/${bam}*.bed*
rm $outdir_endmotif/${bam}*.bam*
rm $outdir_endmotif/${bam}_metrics.txt
rm $outdir_endmotif/${bam}*.bed*
else
echo -e "Errors in end motif run"
fi

echo -e "Script completed running"

EOF

done
