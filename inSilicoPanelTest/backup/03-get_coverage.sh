# define input output files
outDir=$2
bamFile=$1

#basename=$(basename $2)
#fileOut=${basename%.*}
bedFile="$outDir/Panel.bed"




# check that input file exists
if [ -e $bamFile ] && [ ! -z $bamFile ]
then

	# get name for outFile
	fullFile=$(basename "$bamFile")
	fileName=${fullFile%.*}
	outFile="$outDir/$fileName.xls"

	#create out file names
	tmpFile="$outDir/tmp"
	covOut="$outDir/$fileName.cov.xls"
	badCovOut="$outDir/$fileName.badcov.xls"


	#outFile="$outDir/001.xls"

	# create output file
	if [ ! -d $outDir ]
	then
		mkdir -p $outDir
	fi


	#check if bam file has index, else create index
	if [ ! -f $bamFile.bai ]
	then
		echo ""
		echo "#### indexing file ####"
		echo ""
		/data/ngs/bin/samtools/1.3.1/bin/samtools \
		index \
		$bamFile
	fi

	#########
	######## extract coverage according to bedfile locations
	########
	echo ""
	echo "#### extracting coverage $fileName #### "
	/data/ngs/bin/samtools/1.3.1/bin/samtools \
	bedcov \
	--reference /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
	$bedFile \
	$bamFile > $tmpFile
#	$bamFile | tee $tmpFile

	 

	# get relative coverage
	# get file in name
	fileIn=$tmpFile




	# for each line get mean coverage and add to file
	# depth / ( end - start)
	cat $tmpFile | awk '{ meanDepth=($6 / ( $3-$2 ))} {roundMean=sprintf("%.0f", meanDepth)} { print $0 "\t" roundMean}' > $outFile

	# rearrange collumns
	# add header
	echo "#geneName	exonNumber	chr	exonStart	exonEnd	cumCov	meanCov" > $tmpFile
	cat $outFile | awk '{ print $4"\t"$5"\t"$1"\t"$2"\t"$3"\t"$6"\t"$7 }' >> $tmpFile
	mv $tmpFile  $outFile



	# get percent of coverage of panel exones with decent (>= 20 coverage)
	failed=$(cat $outFile | awk '{ if ($7 < 20 ) print}' | wc -l)
	total=$(cat $outFile | wc -l)
	coverage=$(echo "1- $failed/$total" | bc -l)


	#clear screen to make coverage better visible
#	clear
	echo "Coverage = $coverage" | tee $covOut

	head -n 1 $outFile > $badCovOut
	cat $outFile | awk '{ if ($7 < 20) print }' >> $badCovOut

else 
echo "Input \"$bamFile\" file not correct"
fi


### get all nongood exons
echo "#geneName       exonNumber      chr     exonStart       exonEnd" > $outDir/badCoveredExons.xls

for i in $(ls Panel/*badcov*)
do
	tail -n +2 $i
done | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | sort | uniq >> $outDir/badCoveredExons.xls



# clear intermediate file





