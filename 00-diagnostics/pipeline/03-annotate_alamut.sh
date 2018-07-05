#!/bin/bash
#SBATCH --mailType=FAIL
#SBATCH --cpus-per-task=6


#########################
######## prepare ########
#########################

fileIn=$1
db=$2
runAlamut=$3

CPU=6


#### programs
gemini="/data/programs/bin/ngs/gemini/anaconda/bin/gemini"

#fileIn="00-vcf/consensus.vcf"
#db="08-database/consensus-vep-snpEff.db"

#### prepare out directory
alaOut="10-alamut"
mkdir -p $alaOut



###############################################
######## annotate variants with alamut ########
###############################################
bn=$(basename $fileIn)
fn=${bn%.*}
out="$alaOut/$fn.alamut"

## check if alamut needs to be run, or if this has be done before (userdefined)
if [ "$runAlamut" != "false" ]
then
	/data/programs/bin/ngs/alamut/batch/1.8/alamut-batch \
	--outputannonly \
	--in $fileIn \
	--hgmdUser hoppmann \
	--hgmdPasswd KiKli2016 \
	--ann $out \
	--unann $alaOut/unannotated.log
fi


#############################################################
######## prepareing bed file from alamut annotations ########
#############################################################

echo "preparing bed file"
fileIn=$out
bedFile="$alaOut/$fn.bed.gz"

cat $fileIn | awk '{
	if (FNR == 1) {
		gsub("\t", "\talamut_");
		print "#CHROM\tSTART\tSTOP\t" $0
	} else {
		print $2 "\t" ($3 - 1) "\t" $3 + length($4) - 1 "\t" $0
	}
}' | cut -f 1-4,9-999 | sed 's/#id/rsId/' | bgzip > $bedFile

tabix -p bed $bedFile




######## prepare variables to read alamut out in gemini ########

echo "preparing variables"

## extract indexes, coltypes and colmodes
ALAMUT_IDXS=""
ALAMUT_COLTYPES=""
ALAMUT_MODE=""


start=4
nAnno=$(zcat $bedFile | head -n 1 | wc | awk '{ print $2 }')
for (( i=$start; i<=$nAnno; i++ ))
do
	if [ $i -eq $start ]
	then
		ALAMUT_IDXS=$i
		ALAMUT_COLTYPES="text"
		ALAMUT_MODE="list"
	else
		ALAMUT_IDXS="$ALAMUT_IDXS,$i"
		ALAMUT_COLTYPES="$ALAMUT_COLTYPES,text"
		ALAMUT_MODE="$ALAMUT_MODE,list"
	fi
done

#echo $ALAMUT_IDXS
#echo $ALAMUT_COLTYPES
#echo $ALAMUT_MODE

#### get all header besides chrom start stop
ALAMUT_COLS=$(zcat $bedFile | head -n 1 | cut -f $start-999 | tr "\t" ",")
#echo $ALAMUT_COLS
#echo $ALAMUT_IDXS
#echo $ALAMUT_COLTYPES
#echo $ALAMUT_MODE


echo "Updating database"
$gemini annotate $db \
-f $bedFile \
-a extract \
-c $ALAMUT_COLS \
-e $ALAMUT_IDXS \
-t $ALAMUT_COLTYPES \
-o $ALAMUT_MODE






















echo "Script finished"
