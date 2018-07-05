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
runAlamut=$2
# panels=${@:2}



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
/data/programs/bin/ngs/vcftools/0.1.13/bin/vcf-merge $files > $mergedOut




##### annotate variants using VEP and SNPEff
/data/programs/scripts/hoppmann/00-diagnostics/pipeline/01-annotate.sh $mergedDir $annotateDir


##### load annotated files in gemini DB
/data/programs/scripts/hoppmann/00-diagnostics/pipeline/02-load-in-gemini.sh $annotateDir $databaseDir


######## make gemini databases usable
chmod 755 $databaseDir/*.db


###### annotate with Alamut and update database
echo "annotating with alamut"
/data/programs/scripts/hoppmann/00-diagnostics/pipeline/03-annotate_alamut.sh $mergedOut $databaseDir/merged.db $runAlamut


######## make gemini databases usable
chmod 755 $databaseDir/*.db
