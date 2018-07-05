#!/bin/bash

clear

# get dirname of current script directory
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


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
	$scriptDir/02-universal_filter.pl -d $db -pat $patID -pan $i
done


## clean script
$scriptDir/04-clean_results.sh 09-filter/$patID


# add additional inforation about number of predictions
for i in $(ls 09-filter/$patID/*.out)
do

	## exclude files already modyfied
	if [[ $i != *"mod.out" ]]
	then
		in=$i
		out=${in/.out/_mod.out}

		Rscript $scriptDir/05-postFilterResultModifications.R $in $out

		mv $out $in

	fi
done






echo $DIR
