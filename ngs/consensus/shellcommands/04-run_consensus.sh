clear

#~/scripts/ngs/consensus/01-get_consensus.pl \
#-vcfIn test_small.vcf \
#-out bla.vcf



#~/scripts/ngs/consensus/01-get_consensus.pl \
#-vcfIn all.vcf \
#-minHits 1 \
#-out consensus_1.vcf


#~/scripts/ngs/consensus/01-get_consensus.pl \
#-vcfIn all.vcf \
#-minHits 2 \
#-out consensus_2.vcf

#~/scripts/ngs/consensus/01-get_consensus.pl \
#-vcfIn all.vcf \
#-minHits 3 \
#-out consensus_3.vcf

#~/scripts/ngs/consensus/01-get_consensus.pl \
#-vcfIn all.vcf \
#-minHits 4 \
#-out consensus_4.vcf


dirIn="02-merged"
dirOut="03-consensus/oldReadgroups"

mkdir -p $dirOut


for i in $(ls $dirIn/*.vcf)
do
	
	fileName=$(basename $i)
	
	echo "Processing $fileName"
	
	~/scripts/ngs/consensus/01-get_consensus.pl \
	-vcfIn $i \
	-minHits 2 \
	-out $dirOut/$fileName
	
	
done
