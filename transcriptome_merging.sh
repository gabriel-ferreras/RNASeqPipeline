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

PER=$(grep Missed exons: comparison.stats | awk '{ print $2 }')
echo "      Percentage of similiraty of whole transcriptome assembly to reference genome is "$PER

cd ../samples
i=1
if [ $PER -lt 5 ]
then
	while [ $i -le $NUM_SAMPLES ]
	do
        	cd sample_$i
        	stringtie -e -B -G ../../annotation/annotation.gtf -o sample_$i.gtf sample_$i.bam
        	cd ..
        	((i++))
	done
else
	while [ $i -le $NUM_SAMPLES ]
        do
                cd sample_$i
                stringtie -e -B -G ../../results/stringtie_merged.gtf -o sample_$i.gtf sample_$i.bam
                cd ..
                ((i++))
        done
fi

echo "Analysis complete :)"

