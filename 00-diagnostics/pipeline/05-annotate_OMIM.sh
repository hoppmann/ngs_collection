#!/bin/bash
geminiDB=$1
pipelineDir="/data/programs/scripts/hoppmann/00-diagnostics/pipeline/"
outDir="/data/public_resources/geneInfos/BED/"
bedFileName="OMIM.bed"

##########################################
######## prepare and run intervar ########
##########################################

echo "Updating Gemini with Alamut phenotypes"

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


echo "OMIM annotation done!"
