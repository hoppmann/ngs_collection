# as input use a list of genes one gene per row
# recieve file containing genes to check
panelFile=$1
basename=$(basename $panelFile)
fileName=${basename%.*}
outDir=$fileName

mkdir -p $outDir


/h/hoppmann/scripts/ngs/inSilicoPanelTest/01-extract_transcripts.sh $1 $outDir

/h/hoppmann/scripts/ngs/inSilicoPanelTest/02-form_bed.pl $outDir

for i in $(ls 00-bam/*.bam)
do
	. 03-get_coverage.sh $i $outDir
done	
