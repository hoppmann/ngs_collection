#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6

clear

#### recive folder containing vcf files

if [ -z "$1" ]
then
	echo "No folder given. Please specify folder containing vcf-files."
	exit
fi

vcfDir=$1
panels=${@:2}



## define folder variables
mergedDir="01-merged"
annotateDir="02-annotated"
databaseDir="03-databases"
filterDir="04-filter"

mkdir -p $mergedDir


###### merge vcf-files
files=""
for i in $(ls $vcfDir/*.vcf)
do
	echo $i
	bgzip -c $i > "$i.gz"
	tabix "$i.gz"
	files="$files $i.gz"

done


##### merge vcf-files
mergedOut="$mergedDir/merged.vcf"
/data/ngs/bin/vcftools/0.1.13/bin/vcf-merge $files > $mergedOut




##### annotate variants using VEP and SNPEff
~/scripts/ngs/diagnostics/pipeline/01-annotate.sh $mergedDir $annotateDir


##### load annotated files in gemini DB
~/scripts/ngs/diagnostics/pipeline/02-load-in-gemini.sh $annotateDir $databaseDir


######## make gemini databases usable
chmod 755 $databaseDir/*.db


###### annotate with Alamut and update database
echo "annotating with alamut"
~/scripts/ngs/diagnostics/pipeline/03-annotate_alamut.sh $mergedOut $databaseDir/merged.db


######## make gemini databases usable
chmod 755 $databaseDir/*.db



##### filter variants
#for DB in $(ls $databaseDir/*.db)
#do
#	for curPanel in ${panels[@]}
#	do
#		outDir=$(basename $DB)
#		outDir="${outDir%.*}"
#		echo $outDir
#		/h/hoppmann/scripts/ngs/diagnostics/pipeline/03-universal_filter.pl -d $DB -p $curPanel -t
#	done
#done
