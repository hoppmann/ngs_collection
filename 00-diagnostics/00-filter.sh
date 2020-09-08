#!/bin/bash

# clear

## User input variables
patID=$1
shift
# panelDir=$@


# get dirname of current script directory & current filter dir

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pipelineDir="$scriptDir/01-filter/"
resultDir="09-filter"
patResDir="$resultDir/$patID"
panels=$@
panels+=('exome')


if [ -z $patID ]
then
	echo -e "######## ERROR: patID not given\n00-filter.sh $panelDir $patID"
fi


for panelDir; do
	if [ -z $panelDir ]
	then
		echo -e "######## ERROR: panel not given\n00-filter.sh $panelDir $patID"
	fi
done




#################################
######## filter variants ########
#################################

# run for all given panels
for panelDir in ${panels[@]};
do

	## choose correct path for database to use
	db1="08-database/consensus-vep-snpEff.db"
	db2="09-database/consensus-vep-snpEff.db"
	db3="03-databases/merged.db"


	if [ -f $db1 ]; then
		db=$db1
	elif [ -f $db2 ]; then
		db=$db2
	elif [[ -f $db3 ]]; then
		db=$db3
	fi


	if [ "$panelDir" == "exome" ]
	then
		# run in exome mode
		$pipelineDir/02-universal_filter.pl -d $db -pat $patID -outDir $resultDir
	else
		#### run filtering step for panel
		for i in $(ls $panelDir/*);
		do
			$pipelineDir/02-universal_filter.pl -d $db -pat $patID -pan $i -outDir $resultDir
		done
	fi

done





########################################
######## post filter processing ########
########################################


for panelDir in ${panels[@]}
do

	panelName=$(basename $panelDir)

	inDir=$patResDir/$panelName


	# clean script
	$pipelineDir/04-clean_results.sh $inDir



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


			Rscript $pipelineDir/05-postFilterResultModifications.R $in $out

		fi

	done



	### create single excel sheet of all _mod files
	echo "preparing excel file"
	sheets=""
	for i in $(ls $inDir/*_mod.*)
	do

		sheets="$sheets $i"

	done
	bn=$(basename $inDir )


	~/.local/bin/csv2xls \
	$sheets \
	-n \
	-o $inDir/$bn.xls \
	-d '	'



	# ###################################
	# ######## coverage analysis ########
	# ###################################
	#
	# if [ $panelDir != "exome" ]
	# then
	#
	# 	#### check coverage
	# 	for panelFile in $(ls $panelDir/*)
	# 	do
	# 		$pipelineDir/03-coverage_analysis.sh $patID $panelFile $resultDir
	# 	done
	# fi

done

















echo $DIR
