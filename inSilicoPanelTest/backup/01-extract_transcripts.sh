######### extract transcript reagions from database ########




# get file in name
fileIn=$1
file=$(basename $fileIn)
fileName=${file%.*}

# define out dir 
outDir="$2"
fileOut="$outDir/$fileName-transcripts.pos"


# for each gene in gene list extract transcripts
dbIn="00-data/exon.db"
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



# make header
#echo "gene	chr	cds_from	cds_to	cds_locations" > $fileOut

count=0
for i in $(cat $fileIn)
do
	length=$(cat $fileIn | wc -l)
	count=$((count+1))
	echo "$count/$length"
	sqlite3 $dbIn "select gene, chr, cds_from, cds_to, cds_locations from $table where gene=='$i'" >> $fileOut
done 

cat $fileOut | sed 's/|/\t/g'| sed 's/\[//g' | sed 's/]//g' | sed 's/, /\t/g' > $outDir/tmp

mv $outDir/tmp $fileOut

