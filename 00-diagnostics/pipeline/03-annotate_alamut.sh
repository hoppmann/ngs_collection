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


#### check if alamut needs to be run, or if this has be done before (userdefined)
if [ "$runAlamut" != "false" ]
then
	echo "Running Alamut"
	/data/programs/bin/ngs/alamut/batch/1.11/alamut-batch \
	--outputannonly \
	--in $fileIn \
	--ann $out \
	--dogenesplicer \
	--donnsplice \
	--unann $alaOut/unannotated.log

fi

       # --hgmdUser hoppmann \
       # --hgmdPasswd KiKli2016 \


#############################################################
######## prepareing bed file from alamut annotations ########
#############################################################

echo "Preparing bed file"
fileIn=$out
bedFile="$alaOut/$fn.bed.gz"

cat $fileIn | awk '{
	if (FNR == 1) {
		gsub("\t", "\talamut_");
		print "#CHROM\tSTART\tSTOP\t" $0
	} else {
		print $2 "\t" ($3 - 1) "\t" $3 + length($4) - 1 "\t" $0
	}
}' | cut -f 1-3,9-999 | bgzip > $bedFile

# | cut -f 1-4,9-999 | sed 's/#id/rsId/' | bgzip > $bedFile

tabix -p bed $bedFile





###################################################
######## preparing gemini annotate command ########
###################################################



######## prepare variables to read alamut out in gemini ########

echo "Preparing variables"

## extract indexes, coltypes and colmodes
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


start=3


#### get all header besides chrom start stop
ALAMUT_COLS=$(zcat $bedFile | head -n 1 | cut -f $start-999 | tr "\t" ",")
IFS=',' read -r -a ALAMUT_ARRAY <<< "$ALAMUT_COLS"




nAnno=$(zcat $bedFile | head -n 1 | wc | awk '{ print $2 }')

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

#echo $ALAMUT_IDXS
#echo $ALAMUT_COLTYPES
#echo $ALAMUT_MODE


# echo $ALAMUT_COLS
# echo $ALAMUT_IDXS
# echo $ALAMUT_COLTYPES
# echo $ALAMUT_MODE


echo "Updating gemini DB"
$gemini annotate $db \
-f $bedFile \
-a extract \
-c $ALAMUT_COLS \
-e $ALAMUT_IDXS \
-t $ALAMUT_COLTYPES \
-o $ALAMUT_MODE







$gemini annotate \
-f $bedFile \
-c inAlamut \
-a boolean \
$db

































echo "Script finished!"
# done
