#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

INPUTDIR=/cluster/projects/tcge/cell_free_epigenomics/processed_data/TCGE-CFMe-HCC/dedup_bam_umi_pe
base=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/Run5_HCC
outdir_delfi=$base/Output_Delfi
shdir=$base/sh_scripts
fragmentomics_delfi=/cluster/projects/pughlab/bin/fragmentomics/v2/ratio # The DELFI plots

ref=/cluster/projects/tcge/DB/MEDIPIPE/hg38/bwa_index_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
picard_dir=/cluster/tools/software/picard/2.10.9
fragmentomics_peak=/cluster/projects/pughlab/bin/fragmentomics/v2/nucleosome_peak
peaks=/cluster/projects/pughlab/bin/fragmentomics/v2/nucleosome_peak/peaks/hg38
outdir_peak=$base/Output_Peak

fragmentomics_endmotif=/cluster/projects/tcge/cell_free_epigenomics/processed_data/fragmentomics_dory/updated_endmotif_scripts
#fragmentomics_endmotif=/cluster/projects/pughlab/bin/fragmentomics/v2/end_motif/R
outdir_endmotif=$base/Output_Endmotif
outdir_endmotif_all_frags=$base/Output_Endmotif_All_Frags


insertsize=$base/Insert_size

slurm=$base/slurm

mkdir -p $base
mkdir -p $outdir_delfi
mkdir -p $outdir_peak
mkdir -p $outdir_endmotif
mkdir -p $outdir_endmotif_all_frags
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
#SBATCH --time=1-12:00:00
#SBATCH --mem=30GB
#SBATCH -c 1

module load picard/2.10.9
module load igenome-human/hg38
module load samtools
module load R/4.0.0
module load bedtools/2.27.1


#### Run fragmentation ratio

Rscript $fragmentomics_delfi/runFrag.R\
 --id $bam\
 --bamdir $INPUTDIR\
 --filters $fragmentomics_delfi/extdata/filters.hg38.rda\
 --gaps $fragmentomics_delfi/extdata/gaps.hg38.rda\
 --VNTRs $fragmentomics_delfi/extdata/VNTRs.hg38.rda\
 --tiles $fragmentomics_delfi/extdata/hg38_tiles.bed\
 --healthy $fragmentomics_delfi/extdata/healthy.median.hg38.rda\
 --outdir $outdir_delfi\
 --libdir $fragmentomics_delfi
 
#### Run nucleosome peak

module load samtools/1.10
module load picard/2.10.9
module load bedtools/2.27.1




## Subset the bam file for insert sizes of 167bp
samtools view -h $INPUTDIR/${bam}.bam | awk 'substr(\$0,1,1)=="@" || (\$9==167) || (\$9==-167)' | \
samtools view -b > $outdir_peak/${bam}.bam
samtools index $outdir_peak/${bam}.bam

## Remove duplicates
java -Xmx28g -jar $picard_dir/picard.jar MarkDuplicates \
I=$outdir_peak/${bam}.bam \
O=$outdir_peak/${bam}_deduped.bam \
M=$outdir_peak/${bam}_metrics.txt \
TMP_DIR=$outdir_peak \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true \
MAX_RECORDS_IN_RAM=100000 

samtools sort -n $outdir_peak/${bam}_deduped.bam -o $outdir_peak/${bam}_deduped_sorted.bam

## Sort bam by name and convert to bedpe format
samtools view -bf 0x2 $outdir_peak/${bam}_deduped_sorted.bam | \
bedtools bamtobed \
-i stdin \
-bedpe > $outdir_peak/${bam}.bedpe

## Split bedpe into chromosomes
$fragmentomics_peak/R/splitBedpe.sh $outdir_peak/${bam}.bedpe

## Calculate the distance from closest peak
Rscript $fragmentomics_peak/R/nucleosome_peaks_distance.R \
--id $bam \
--path $outdir_peak \
--peaks $peaks \
--outdir $outdir_peak

#### Run end motif - continuing on last steps, 167bp fragments

## Format bedpe to bed
Rscript $fragmentomics_endmotif/motif_format_bedpe.R \
--id $bam \
--bedpe $outdir_peak/${bam}.bedpe \
--outdir $outdir_endmotif

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

## Get insert size 
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
INPUT=$INPUTDIR/${bam}.bam OUTPUT=$insertsize/${bam}.insert_size.txt H=$insertsize/${bam}.insert_size_histogram.pdf M=0.05


#### Cleaning up
## Remove intermediate files - nucleosome peak
if [[ -f "$outdir_peak/${bam}_peak_distance.txt" ]]
then
echo -e "Nucleosome peak output completed sucessfully"
rm $outdir_peak/${bam}*.bam*
rm $outdir_peak/${bam}_metrics.txt
rm $outdir_peak/${bam}*.bed*
else
echo -e "Errors in nucleosome peak run"
fi

## Remove intermediate files - end motif
if [[ -f "$outdir_endmotif/${bam}_raw.txt" ]]
then
echo -e "Endmotif output completed sucessfully on all fragments 20 to 600bp"
rm $outdir_endmotif/${bam}*.bam*
rm $outdir_endmotif/${bam}_metrics.txt
rm $outdir_endmotif/${bam}*.bed*
else
echo -e "Errors in end motif run on 167bp fragments"
fi


#### Now rerun end motifs on all fragments

## Subset the bam file for insert sizes between 20bp and 600bp
samtools view -h $INPUTDIR/${bam}.bam | awk 'substr(\$0,1,1)=="@" || (\$9>=20 && \$9<=600) || (\$9<=-20 && \$9>=-600)' | \
samtools view -b > $outdir_peak/${bam}.bam
samtools index $outdir_peak/${bam}.bam

## Remove duplicates
java -Xmx28g -jar $picard_dir/picard.jar MarkDuplicates \
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

#### Run end motif - continuing on last steps all fragments
## Format bedpe to bed
Rscript $fragmentomics_endmotif/motif_format_bedpe.R \
--id $bam \
--bedpe $outdir_peak/${bam}.bedpe \
--outdir $outdir_endmotif_all_frags

rm $outdir_peak/${bam}.bedpe

## Get FASTA sequences 
bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir_endmotif_all_frags/${bam}_5.bed > $outdir_endmotif_all_frags/${bam}_fasta_5.bed
bedtools getfasta \
-bedOut \
-fi $ref \
-bed $outdir_endmotif_all_frags/${bam}_3.bed > $outdir_endmotif_all_frags/${bam}_fasta_3.bed

## Convert FASTA to end motif context frequencies
Rscript $fragmentomics_endmotif/motif_get_contexts.R \
--id $bam \
--fasta_5 $outdir_endmotif_all_frags/${bam}_fasta_5.bed \
--fasta_3 $outdir_endmotif_all_frags/${bam}_fasta_3.bed \
--outdir $outdir_endmotif_all_frags


## Remove intermediate files - end motif
if [[ -f "$outdir_endmotif_all_frags/${bam}_raw.txt" ]]
then
echo -e "Endmotif output completed sucessfully on all fragments 20 to 600bp"
rm $outdir_endmotif_all_frags/${bam}*.bam*
rm $outdir_endmotif_all_frags/${bam}_metrics.txt
rm $outdir_endmotif_all_frags/${bam}*.bed*
else
echo -e "Errors in end motif run on 20 to 600bp fragments"
fi


echo -e "Script completed running"

EOF

done
