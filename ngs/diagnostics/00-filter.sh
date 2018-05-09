#!/bin/bash

clear

db="08-database/consensus-vep-snpEff.db"
#db="03-databases/merged.db"
panelDir=$1
patID=$2

if [ -z $panelDir ]
then
	echo "######## ERROR: panel not given"
fi

if [ -z $patID ]
then
	echo "######## ERROR: patID not given"
fi


for i in $(ls $panelDir/*);
do
	~/scripts/ngs/diagnostics/02-universal_filter.pl -d $db -pat $patID -pan $i
done


~/scripts/ngs/diagnostics/04-clean_results.sh 09-filter/$patID
