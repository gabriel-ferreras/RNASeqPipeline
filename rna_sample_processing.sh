#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

SAMPLE_DIR=$1
i=$2
NUM_SAMPLES=$3

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome
fastqc sample_$i.fastq.gz
hisat2 --dta -x ../../genome/index -U sample_$i.fastq.gz -S sample_$i.sam

## Generating sorted bam file
samtools sort -o sample_$i.bam sample_$i.sam
rm sample_$i.sam
rm *.fastq
samtools index sample_$i.bam

## Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample_$i.gtf -l sample_$i sample_$i.bam

## Preparing merge list file for transcriptome merging.
echo $SAMPLE_DIR/sample_$i.gtf > ../../results/merge_list.txt

## Communication with blackboard.
echo "Processing Sample $i done!" >> ../../results/blackboard.txt
NUM_PROC=$(wc -l ../../results/blackboard.txt)

if [ $NUM_PROC -eq $NUM_SAMPLES ]
do

done
