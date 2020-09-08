
##################################
######## gather variables ########
##################################
inputFile=$1
outDir=$2
geminiDB=$3
pipelineDir="/data/programs/scripts/hoppmann/00-diagnostics/pipeline/"
intervarConfigFile="interVar_config.ini"
bedFileName="interVar.bed"

outFileName=$(basename $inputFile)
outFileName=${outFileName%%.vcf}
outFile="$outDir/$outFileName"


##########################################
######## prepare and run intervar ########
##########################################

mkdir -p $outDir
cp $pipelineDir/$intervarConfigFile $outDir/$intervarConfigFile

sed -i "s#{INPUTFILE}#$inputFile#" $outDir/$intervarConfigFile
sed -i "s#{OUTPUTFILE}#${outFile}#" $outDir/$intervarConfigFile







python /data/programs/bin/ngs/InterVar/Intervar.py  -c $outDir/$intervarConfigFile




###########################################
######## annotate Gemini with ACMG ########
###########################################


### prepare BED file
echo "Prepareing BED"
cat $outFile.hg19_multianno.txt.intervar | awk '{
	if (FNR == 1) {
		gsub("\t", "\t");
		print "#CHROM\tSTART\tSTOP\t" $0
	} else {
		print $1 "\t" ($3 - 1) "\t" $3 + length($4) - 1 "\t" $0
	}
}' | cut -f 1-3,17 | sed 's/PVS1/\tPVS1/' | sed 's/ and /\t/' | sed 's/InterVar: //' | sed 's/InterVar/ACMG/' > $outDir/$bedFileName


#### sort bed file for tabix usage
echo "Sort BED-file"
/data/programs/bin/ngs/bedtools/2.26.0/bin/sortBed -header -i $outDir/$bedFileName | bgzip > $outDir/$bedFileName.gz


#### index .bed.gz file
tabix -f -p bed $outDir/$bedFileName.gz


########################################
######## add to gemini database ########
########################################


echo "Updating Gemini"

ALAMUT_IDXS=""
ALAMUT_COLTYPES=""
ALAMUT_MODE=""



function chooseType() {

	if [[ $1 == *gnomad* ]]; then
		colMode="max"
		colType="float"
	elif [[ $1 == *1000g* ]]; then
		colMode="max"
		colType="float"

	elif [[ $1 ==  *esp* ]]; then
		colMode="max"
		colType="float"

	elif [[ $1 == *rsId* ]]; then
		colMode="first"
		colType=text
	else
		colMode="list"
		colType="text"
	fi

}


start=4



#### get all header besides chrom start stop
ALAMUT_COLS=$(zcat $outDir/$bedFileName.gz | head -n 1 | cut -f $start-999 | tr "\t" ",")
IFS='\t' read -r -a ALAMUT_ARRAY <<< "$ALAMUT_COLS"




nAnno=$(zcat $outDir/$bedFileName.gz | head -n 1 | tr " " "_" | wc | awk '{ print $2 }')

for (( i=$start; i<=$nAnno; i++ ))
do

	chooseType "${ALAMUT_ARRAY[$i - 4]}"

	if [ $i -eq $start ]
	then
		ALAMUT_IDXS=$i
		ALAMUT_COLTYPES=$colType
		ALAMUT_MODE=$colMode
	else
		ALAMUT_IDXS="$ALAMUT_IDXS,$i"
		ALAMUT_COLTYPES="$ALAMUT_COLTYPES,$colType"
		ALAMUT_MODE="$ALAMUT_MODE,$colMode"
	fi
done

gemini annotate  \
-f $outDir/$bedFileName.gz \
-a extract \
-c $ALAMUT_COLS \
-e $ALAMUT_IDXS \
-t $ALAMUT_COLTYPES \
-o $ALAMUT_MODE \
$geminiDB


echo "ACMG annotation done!"






















































#
