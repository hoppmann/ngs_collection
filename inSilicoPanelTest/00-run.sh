# as input use a list of genes one gene per row
# recieve file containing genes to check
panelFile=$1
bamDir=$2

## check if panel file is given and exists

if [ -z $panelFile ] || [ ! -f $panelFile ]
then
	echo "Panel file not given or doesn't exist."
	exit
fi

basename=$(basename $panelFile)
fileName=${basename%.*}
outDir=$fileName

#bamDir="00-bam"



# check if outDir exists. if so delete it
if [ -d $outDir ]
then
	rm -r $outDir
fi


maxCPU=$(cat /proc/cpuinfo | grep processor  | wc -l)



mkdir -p $outDir

/h/hoppmann/scripts/ngs/inSilicoPanelTest/01-extract_transcripts.sh $panelFile $outDir


echo "Creating bed file"
/h/hoppmann/scripts/ngs/inSilicoPanelTest/02-form_bed.pl $outDir $fileName



for i in $(ls $outDir/bed/*.bed)
do

	# init flag variable
	indexFiles=0

	# check if there is an indexed bam file
	for j in $(ls $bamDir/*.bam)
	do
		if [ ! -e $j.bai ]
		then
			#if an index is missing make flag
			indexFiles=1
		fi
	done

	# if there are unindexed files run single bam else multithred
	if [ $indexFiles == 1 ]
	then
		echo "Get coverage: $i non multithred"
		~/scripts/ngs/inSilicoPanelTest/03-coverage.pl $i $bamDir $outDir
	else
		echo "Get coverage: $i; multithreaded"
		~/scripts/ngs/inSilicoPanelTest/03-coverage.pl $i $bamDir $outDir &
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
done > $outDir/00-CoverageSummary.xls
