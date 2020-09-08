# as input use a list of genes one gene per row
# recieve file containing genes to check
bamFile=$1
panelFile=$2
resultDir=$3
baiFile=${bamFile/.bam/.bai}



outName="00-BadCoveredExons.txt"
# resultDir="09-filter"

## add path to pipeline
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
pipelinePath=`dirname $SCRIPT`

## check if panel file is given and exists

if [ -z $panelFile ] || [ ! -f $panelFile ]
then
	echo "Panel file not given or doesn't exist."
	exit
fi

basename=$(basename $panelFile)
fileName=${basename%.*}

basename=$(basename $bamFile)
patID=${basename%.*}


covOutDir="$resultDir/$patID/00-coverage/"
outDir="$covOutDir/$fileName"



# check if outDir exists. if so delete it
if [ -d $outDir ]
then
	rm -r $outDir
fi


maxCPU=$(cat /proc/cpuinfo | grep processor  | wc -l)



mkdir -p $outDir

# extract transcripts
$pipelinePath/01-extract_transcripts.sh $panelFile $outDir


echo "Creating bed file"
$pipelinePath/02-form_bed.pl $outDir $fileName



for i in $(ls $outDir/bed/*.bed)
do

	# init flag variable
	indexFiles=0

	# check if there is an indexed bam file
	if [ ! -e $baiFile ]
	then
		#if an index is missing make flag
		indexFiles=1
	fi

	# if there are unindexed files run single bam else multithred
	if [ $indexFiles == 1 ]
	then
		echo "Get coverage: $i"
		$pipelinePath/03-coverage_singlebam.pl $i $bamFile $outDir
	else
		echo "Get coverage: $i"
		$pipelinePath/03-coverage_singlebam.pl $i $bamFile $outDir &
	fi


	# make scatteler
	while [ $(jobs | wc -l) -gt $maxCPU ]
	do
		sleep 0.01
	done


done

wait




# extract overall coverage for each out file
if [ -e $outDir/$outName ]
then
	rm $outDir/$outName
fi

for i in $(ls $outDir/*.txt)
do
	# echo $i
	if [ $(awk '/Bad covered exons/,/^$/' $i | tail -n +3 | sed '/^$/d'| wc -l) -gt 0 ]
	then
		fn=$(basename -- $i)
		fn=${fn%.*}
		echo $fn

		head -n 1 $i
		awk '/Bad covered exons/,/^$/' $i
		echo ""
	fi


done > $covOutDir/$fileName.txt

rm -r $outDir



















































# echo "Done"
