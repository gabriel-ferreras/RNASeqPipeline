#! /bin/bash

#Print help message with usage if no parameters indicated.
if [ $# -ne 1 ]
then
    echo ""
    echo "    Usage: bash RNASeqPipeline.sh <params_file>"
    echo ""
    echo "    params_file: imput file with the parameters."
    echo ""
    exit
fi

#Reading parameters file.
PARAMS=$1
echo ""
echo "======================"
echo "| LOADING PARAMETERS |"
echo "======================"
echo ""
WORK_DIR=$(grep working_directory: $PARAMS | awk '{ print $2 }')
echo "      Working directory = "$WORK_DIR

EXP=$(grep experiment_name: $PARAMS | awk '{ print $2 }')
echo "      Experiment name = "$EXP

NUM_SAMPLES=$(grep number_samples: $PARAMS | awk '{ print $2 }')
echo "      Number samples = "$NUM_SAMPLES

ANNOTATION=$(grep path_annotation: $PARAMS | awk '{ print $2 }')
echo "      Annotation in "$ANNOTATION

GENOME=$(grep path_genome: $PARAMS | awk '{ print $2 }')
echo "      Genome in "$GENOME

SAMPLES=()
i=0
while [ $i -lt $NUM_SAMPLES ]
do
	j=$(( $i + 1 ))
       	SAMPLES[$i]=$(grep path_sample_$j: $PARAMS | awk '{ print $2 }')
	echo "      Sample $j in ${SAMPLES[$i]}"
        ((i++))
done
#Preparing working workspace.
echo ""
echo "======================"
echo "| CREATING WORKSPACE |"
echo "======================"
echo ""
cd $WORK_DIR
mkdir $EXP
cd $EXP
mkdir results genome annotation samples
cd samples
i=1
while [ $i -le $NUM_SAMPLES ]
do
	mkdir sample_$i
	((i++))
done
cd ..

#Copying the data.
cp $ANNOTATION $WORK_DIR/$EXP/annotation/annotation.gtf
cp $GENOME $WORK_DIR/$EXP/genome/genome.fa
i=1
while [ $i -le $NUM_SAMPLES ]
do
	j=$((i - 1))
        cp ${SAMPLES[j]} $WORK_DIR/$EXP/samples/sample_$i/sample_$i.fastqc.gz
        ((i++))
done


#Creating genome index.
echo ""
echo "========================="
echo "| CREATING GENOME INDEX |"
echo "========================="
echo ""
cd $WORK_DIR/$EXP/annotation
extract_splice_sites.py annotation.gtf > splice_sites.ss
extract_exons.py annotation.gtf > exons.exon
cd ../genome
hisat2-build --ss ../annotation/splice_sites.ss --exon ../annotation/exons.exon genome.fa index

#Quality analysis of samples fastq files by FASTQC.
echo ""
echo "=============================="
echo "| QUALITY CONTROL OF SAMPLES |"
echo "=============================="
echo ""

cd ../samples
while [ $i -le $NUM_SAMPLES ]
do
        cd sample_$i
	gzip -d sample_$i.fastqc.gz
	fastqc sample_$i.fastqc
	cd ..
	((i++))
done

#Mapping of lectures to reference genome by HISAT2.
echo ""
echo "========================="
echo "|   MAPPING BY HISAT2   |"
echo "========================="
echo ""



while [ $i -le $NUM_SAMPLES ]
do
	cd sample_$i
	if [ -f SRR5602506_2.fastq ]
then
   fastqc SRR5602506_1.fastq
   fastqc SRR5602506_2.fastq

   hisat2 --dta -x ../../genome/index -1 SRR5602506_1.fastq -2 SRR5602506_2.fastq -S lux6_2.sam
else
   fastqc SRR5602506_1.fastq

   hisat2 --dta -x ../../genome/index -U SRR5602506_1.fastq -S lux6_2.sam
fi
	hisat2 --dta -x ../../genome/index -U sample_$i.fastqc -S sample_$i.sam
	cd ..
 	((i++))
done


echo "Analysis complete :)"
