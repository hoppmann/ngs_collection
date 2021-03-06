#!/bin/bash
#SBATCH --mailType=FAIL
#SBATCH --cpus-per-task=6



## recive in and out directories
inDir=$1
outDir=$2

mkdir -p $outDir

echo "load in gemini"

#define path
# gemini="/data/programs/bin/ngs/gemini/bin/gemini"
gemini="/data/programs/bin/ngs/gemini/anaconda/bin/gemini"
#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)
CPU=20

## remove old db if exists to avoid locked DB error


# old version didn't annotate cdna position due to snpEff
for i in $(ls $1/*.vcf.gz)
do


	filename=$(basename "$i")
	# filename="${filename%-vep.vcf.gz}"
	filename="${filename%-vep-snpeff.vcf.gz}"


	if [ -e $outDir/$filename.db ]
	then
		rm $outDir/$filename.db
	fi



	echo $filename
	echo $i

	$gemini load \
	-t all \
	--cores $CPU \
	-v $i \
	$outDir/$filename.db
done
