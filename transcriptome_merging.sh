#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

RES_DIR=$1
NUM_SAMPLES=$2

## Accessing results folder
cd $RES_DIR

#Whole transcriptome assembly and comparison to reference genome.
echo ""
echo "===================================="
echo "|   WHOLE TRANSCRIPTOME ASSEMBLY   |"
echo "===================================="
echo ""

stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt

echo ""
echo "======================================"
echo "|   COMPARISON TO REFERENCE GENOME   |"
echo "======================================"
echo ""

gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf

#Gene Expression Quantification in each and every sample.
echo ""
echo "===================================="
echo "|  GENE EXPRESSION QUANTIFICATION  |"
echo "===================================="
echo ""

cd ../samples
i=1
while [ $i -le $NUM_SAMPLES ]
do
        cd sample_$i
        #stringtie -e -B -G ../../annotation/annotation.gtf -o sample_$i.gtf sample_$i.bam
        cd ..
        ((i++))
done

echo "Analysis complete :)"

