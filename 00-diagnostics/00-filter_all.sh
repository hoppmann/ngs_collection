#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL

# clear

## User input variables
if [[ ! $1 ]]
then
	resultDir="09-filter"
else
	resultDir="$1"
fi
# shift

#### get dirname of current script directory & current filter dir

# scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
scriptDir="/data/programs/scripts/hoppmann/00-diagnostics"
pipelineDir="$scriptDir/01-filter/"
columnsFile="$pipelineDir/01-columns.sh"
curDir=$(basename $(pwd))
curDir=${curDir##/}

# patResDir="$resultDir/$patID"
# panels=$@
# panels+=('exome')

echo "$curDir"

outDir=$(dirname $resultDir)

#### choose correct path for database to use
if [ $outDir ]; then
	db1="$outDir/08-database/consensus-vep-snpEff.db"
	db2="$outDir/09-database/consensus-vep-snpEff.db"
	db3="$outDir/03-databases/merged.db"
else
	db1="08-database/consensus-vep-snpEff.db"
	db2="09-database/consensus-vep-snpEff.db"
	db3="03-databases/merged.db"
fi

if [ -f $db1 ]; then
	db=$db1
elif [ -f $db2 ]; then
	db=$db2
elif [[ -f $db3 ]]; then
	db=$db3
fi



# get colums
columns=$(cat $columnsFile | grep -v "^#" \
| grep -v "^$" | tr "\n" ", " | sed -r 's/(.*),/\1\n/')



#### filter small set
echo "filtering small set"
maf="0.05"
gemini query -q "select $columns from variants where \
(alamut_gnomadAltFreq_all <= $maf or alamut_gnomadAltFreq_all is null) \
and (alamut_gnomadAltFreq_afr <= $maf or alamut_gnomadAltFreq_afr is null) \
and (alamut_gnomadAltFreq_eas <= $maf or alamut_gnomadAltFreq_eas is null) \
and (alamut_gnomadAltFreq_nfe <= $maf or alamut_gnomadAltFreq_nfe is null) \
and (alamut_gnomadAltFreq_sas <= $maf or alamut_gnomadAltFreq_sas is null) \
and inAlamut == 1 \
and (impact_severity == 'HIGH' or impact_severity == 'MED' or impact_severity is null) \
" \
--header \
$db > $resultDir/$curDir\_small.out





#### filter big set
echo "filtering big set"
gemini query -q " select $columns from variants where \
(alamut_gnomadAltFreq_all <= $maf or alamut_gnomadAltFreq_all is null) \
and (alamut_gnomadAltFreq_afr <= $maf or alamut_gnomadAltFreq_afr is null) \
and (alamut_gnomadAltFreq_eas <= $maf or alamut_gnomadAltFreq_eas is null) \
and (alamut_gnomadAltFreq_nfe <= $maf or alamut_gnomadAltFreq_nfe is null) \
and (alamut_gnomadAltFreq_sas <= $maf or alamut_gnomadAltFreq_sas is null) \
and inAlamut == 1 \
and (impact_severity == 'HIGH' or impact_severity == 'MED' or impact_severity == 'LOW' or impact_severity is null) \
" \
--header \
$db > $resultDir/$curDir\_big.out



#### add additional information about number of predictions

for i in $(ls $resultDir/*.out)
do

	## exclude files already modyfied
	if [[ $i != *"mod.out" ]]
	then
		in=$i
		out=${in/.out/_mod.out}

				## replace ' to avoid errors occuring in R
				sed -i "s/'/_/g" $in

		echo "modifying $in"
		# echo "java -jar $pipelineDir/PostFilterModifications.jar $in $out"
		java -jar $pipelineDir/PostFilterModifications.jar $in $out
		# Rscript $pipelineDir/05-postFilterResultModifications.R $in $out

	fi

done
