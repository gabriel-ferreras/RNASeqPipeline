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

WORK_DIR=$(grep working_directory: $PARAMS | awk '{ print $2 }')
echo "      Working directory = "$WORK_DIR

EXP=$(grep experiment_name: $PARAMS | awk '{ print $2 }')
echo "      Experiment name = "$EXP

NUM_SAMPLES=$(grep number_samples: $PARAMS | awk '{ print $2 }')
echo "      Number samples = "$NUM_SAMPLES

ANNOTATION=$(grep annotation: $PARAMS | awk '{ print $2 }')
echo "      Annotation in "$ANNOTATION

GENOME=$(grep genome: $PARAMS | awk '{ print $2 }')
echo "      Genome in "$GENOME

i=1
while [ $i -le $NUM_SAMPLES ]
do
       	SAMPLE_$i=$(grep sample_$i: $PARAMS | awk '{ print $2 }')
	echo "      Sample $i in = "${SAMPLE_$i}
        ((i++))
done


#Preparing working workspace.
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
mv $ANNOTATION $WORK_DIR/$EXP/annotation/annotation.gtf
mv $GENOME $WORK_DIR/$EXP/genome/genome.fa
i=1
while [ $i -le $NUM_SAMPLES ]
do
        ${SAMPLE_$i} $WORK_DIR/samples/sample_$i/sample_$i.fastqc
        ((i++))
done


#Creating genome index
#cd $WORK_DIR/$EXP/annotation
#extract_splice_sites.py annotation.gtf>splice_sites.ss
#extract_exons.py annotation.gtf>exons.exon
#cd ../genome
#hisat2-build --ss ../annotation/splice_sites.ss --exon ../annotation/exons.exon genoma.fa index



echo "Analysis complete :)"
