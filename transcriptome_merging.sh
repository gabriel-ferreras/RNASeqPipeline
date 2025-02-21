#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

RES_DIR=$1
INS_DIR=$2
EXP_DESIGN=$3

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

## Warning message if the percentage of novel loci is higher than 5%.

# Option 1: with cut.
PER=$(grep "Novel loci:" comparison.stats | awk '{ print $5 }' | cut -d"%" -f 1 | cut -d"." -f 1)

# Option 2: with sed.
#PER=$(grep "Novel loci:" comparison.stats | awk '{ print $5 }' | sed "s/%)//g")

if [ $PER -ge 5 ]
then
	echo "------------------------------------------------------------------------------------------------"
	echo "|WARNING: Percentage of difference of whole transcriptome assembly to reference genome is $PER.|"
	echo "------------------------------------------------------------------------------------------------"
fi

echo "Analysis complete :)"

