#/data/ngs/bin/bedtools/2.26.0/bin/coverageBed \
#-abam IonXpress_001_rawlib.bam \
#-b ../02-extract_exons/panel.bed \
#-counts



#/dsk/localall/lib/java/jre1.8.0_92_x64/bin/java \
#-jar /data/ngs/bin/GATK/3.6.0/GenomeAnalysisTK.jar \
#-T DepthOfCoverage \
#-o outFile \
#-I IonXpress_001_rawlib.bam \
#-R /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
#-geneList ../01-create-gene-panel/Panel_genes.final.csv


# define input output files
outDir="04-coverage"
#bedFile=$1
bedFile="03-extract_exons/panel.bed"
bamFile="00-bam/IonXpress_001_rawlib.bam"
#bamFile="00-bam/IonXpress_002_rawlib.bam"
#bamFile="00-bam/IonXpress_003_rawlib.bam"
#bamFile="00-bam/IonXpress_004_rawlib.bam"

outFile="$outDir/001.xls"

# create output file
if [ ! -d $outDir ]
then
	mkdir -p $outDir
fi


#check if bam file has index, else create index
if [ ! -f $bamFile.bai ]
then
	echo "indexing file"
	/data/ngs/bin/samtools/1.3.1/bin/samtools \
	index \
	$bamFile
fi


# extract coverage according to bedfile locations
echo "extracting coverage"
/data/ngs/bin/samtools/1.3.1/bin/samtools \
bedcov \
--reference /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
$bedFile \
$bamFile | tee $outFile



clear 

# get relative coverage
# get file in name
fileIn=$outFile

# get name for intermediate file
fullFile=$(basename "$fileIn")
fileName=${fullFile%.*}
intermedFile="$outDir/$fileName.mean.xls"


# for each line get mean coverage and add to file
# depth / ( end - start)
cat $fileIn | awk '{ meanDepth=($5 / ( $3-$2 ))} {roundMean=sprintf("%.0f", meanDepth)} { print $0 "\t" roundMean}' > $intermedFile

# get percent of coverage of panel exones with decent (>= 20 coverage)
failed=$(cat $intermedFile | awk '{ if ($6 < 20 ) print}' | wc -l)
total=$(cat $intermedFile | wc -l)
coverage=$(echo "1- $failed/$total" | bc -l)

echo "Coverage = $coverage"
