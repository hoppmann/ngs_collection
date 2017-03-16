#!/bin/bash
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
annotateDir="01-annotated"
databaseDir="02-databases"
filterDir="03-filter"



##### annotate variants using VEP and SNPEff
#~/scripts/ngs/diagnostics/pipeline/01-annotate.sh $vcfDir $annotateDir

##### load annotated files in gemini DB
#~/scripts/ngs/diagnostics/pipeline/02-load-in-gemini.sh $annotateDir $databaseDir


#### make gemini databases usable
chmod 755 $databaseDir/*.db

## filter variants
for DB in $(ls $databaseDir/*.db)
do
	for curPanel in ${panels[@]}
	do
		outDir=$(basename $DB)
		outDir="${outDir%.*}"
		echo $outDir
		/h/hoppmann/scripts/ngs/diagnostics/pipeline/00-universal_filter.pl -d $DB -p $curPanel -t -autoOut
	done
done
