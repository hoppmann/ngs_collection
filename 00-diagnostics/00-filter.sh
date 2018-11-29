#!/bin/bash

clear

## User input variables
patID=$1
shift
# panelDir=$@



# get dirname of current script directory & current filter dir
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
inDir="09-filter/$patID/"


# run for all given panels
for panelDir;
do


	if [ -z $panelDir ]
	then
		echo "######## ERROR: panel not given\n00-filter.sh panelDir patID"
	fi

	if [ -z $patID ]
	then
		echo "######## ERROR: patID not given\n00-filter.sh panelDir patID"
	fi




	## choose correct path for database to use
	db1="08-database/consensus-vep-snpEff.db"
	db2="03-databases/merged.db"

	if [ -f $db1 ]
	then
		db=$db1
	elif [ -f $db2 ]
	then
		db=$db2
	fi





	if [ "$panelDir" == "exome" ]
	then
		# run in exome mode
		$scriptDir/02-universal_filter.pl -d $db -pat $patID
	else
		#### run filtering step for panel
		for i in $(ls $panelDir/*);
		do
			$scriptDir/02-universal_filter.pl -d $db -pat $patID -pan $i
		done
	fi

done




## clean script
$scriptDir/04-clean_results.sh 09-filter/$patID


# add additional information about number of predictions
for i in $(ls $inDir/*.out)
do

	## exclude files already modyfied
	if [[ $i != *"mod.out" ]]
	then
		in=$i
		out=${in/.out/_mod.out}

                ## replace ' to avoid errors occuring in R
                sed -i "s/'/_/g" $in


		Rscript $scriptDir/05-postFilterResultModifications.R $in $out

	fi
done



#### create single excel sheet of all _mod files
echo "preparing excel file"
sheets=""
for i in $(ls $inDir/*_mod.*)
do

	sheets="$sheets $i"

done
bn=$(basename $inDir )


~/.local/bin/csv2xls \
$sheets \
-o $inDir/$bn.xlsx \
-d '	'






echo $DIR
