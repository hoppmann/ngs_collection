#!/bin/bash
#SBATCH --mailType=FAIL
#SBATCH --cpus-per-task=6



## recive in and out directories
inDir=$1
outDir=$2

mkdir -p $outDir



#define path
gemini="/data/programs/bin/ngs/bin/gemini/anaconda/bin/gemini"

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)
CPU=6



# old version didn't annotate cdna position due to snpEff
for i in $(ls $1/*.vcf.gz)
do

        filename=$(basename "$i")
#        filename="${filename%.*}"
        filename="${filename%-vep-snpeff.vcf.gz}"

	echo $filename
	echo $i

	$gemini load \
	-t all \
	--cores $CPU \
  --save-info-string \
	-v $i \
	$outDir/$filename.db
done
