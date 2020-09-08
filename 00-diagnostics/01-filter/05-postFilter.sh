#!/bin/bash

# retrieve variables
inDir=$1

# get dirname of current script directory
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


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
