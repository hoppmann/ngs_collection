######### extract transcript reagions from database ########




# get file in name
fileIn=$1
file=$(basename $fileIn)
fileName=${file%.*}

# define out dir
outDir="$2"
fileOut="$outDir/$fileName-transcripts.pos"


# for each gene in gene list extract transcripts
#scriptPath=$(dirname $0)
#dbIn="$scriptPath/00-data/exon.db"
## add path to pipeline
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
pipelinePath=`dirname $SCRIPT`


dbIn="$pipelinePath/exon.db"
table=CCDS_GRCh37

# create directory if neede
if [ ! -d $outDir ]
then
	mkdir -p $outDir
fi

# remove old file if exists
if [ -f $fileOut ]
then
	rm $fileOut
fi



count=0
for i in $(cat $fileIn | awk '{print $1}')
do
	length=$(cat $fileIn | wc -l)
	count=$((count+1))
	echo "$count/$length"
	sqlite3 $dbIn "select gene, chr, cds_from, cds_to, cds_locations from $table where gene=='$i'" >> $fileOut
done

cat $fileOut | sed 's/|/\t/g'| sed 's/\[//g' | sed 's/]//g' | sed 's/, /\t/g' > $outDir/tmp

mv $outDir/tmp $fileOut
