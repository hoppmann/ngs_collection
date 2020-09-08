# as input use a list of genes one gene per row
# recieve file containing genes to check
bamFile=$1
panelFile=$2

## add path to pipeline
pipelinePath="/data/programs/scripts/hoppmann/00-diagnostics/panelCoverage/"
maxCPU=$(cat /proc/cpuinfo | grep processor  | wc -l)


## check if panel file is given and exists

if [ -z $panelFile ] || [ ! -f $panelFile ]
then
	echo "Panel file not given or doesn't exist."
	exit
fi

basename=$(basename $panelFile)
fileName=${basename%.*}
outDir=$fileName



# check if outDir exists. if so delete it
if [ -d $outDir ]
then
	rm -r $outDir
fi





mkdir -p $outDir

$pipelinePath/01-extract_transcripts.sh $panelFile $outDir


echo "Creating bed file"
$pipelinePath/02-form_bed.pl $outDir $fileName



for i in $(ls $outDir/bed/*.bed)
do

	# init flag variable
	indexFiles=0

	# check if there is an indexed bam file
	if [ ! -e $bamFile.bai ]
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
for i in $(ls $outDir/*.csv)
do
	echo $i
	head -n 1 $i
	echo ""
done > $outDir/00-CoverageSummary.xls
