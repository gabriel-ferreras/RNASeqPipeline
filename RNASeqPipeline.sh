#! /bin/bash

#Print help message with usage if no parameters indicated.
if [ $# -ne 1 ]
then
    echo ""
    echo "    Usage: bash RNASeqPipeline.sh <params_file>"
    echo ""
    echo "    params_file: imput file with the parameters."
    echo "    Example of params_file in test_params."
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


#Processing samples.
echo ""
echo "======================="
echo "|  SAMPLE PROCESSING  |"
echo "======================="
echo ""

cd ../results
touch blackboard.txt
i=1
while [ $i -le $NUM_SAMPLES ]
do
        qsub -o sample_$i -N sample_$i rna_sample_processing $WORK_DIR/$EXP/samples/sample_$i $i $NUM_SAMPLES
        ((i++))
done

#Whole transcriptome assembly and comparison to reference genome.
echo ""
echo "===================================="
echo "|   WHOLE TRANSCRIPTOME ASSEMBLY   |"
echo "===================================="
echo ""

cd ../results
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt
gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf

#Gene Expression Quantification in each and every sample.
echo ""
echo "===================================="
echo "|  GENE EXPRESSION QUANTIFICATION  |"
echo "===================================="
echo ""


i=1
while [ $i -le $NUM_SAMPLES ]
do
        cd sample_$i
        stringtie -e -B -G ../../annotation/annotation.gtf -o sample_$i.gtf sample_$i.bam
        cd ..
        ((i++))
done

echo "Analysis complete :)"
