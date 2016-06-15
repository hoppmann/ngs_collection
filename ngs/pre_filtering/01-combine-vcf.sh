#### merge SNP and INDEL files 
mkdir -p 01-merged

## set path to GATK
GATK="/data/ngs/bin/GATK/3.2-2"


for i in $(ls *snp*.vcf)
 
do
	SNP=$(basename "$i")
	INDEL=${SNP/snp/indel}
	OUT="${INDEL%.*}"
	OUT=${OUT/-indel/}
	echo "Merging $SNP and $INDEL"
	echo "Writing in $OUT"
	
		 
	java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-K $GATK/yong.li_uniklinik-freiburg.de.key \
	-et NO_ET -R /data/ngs/programs/bundle_2.8/ucsc.hg19.fasta \
	--variant $SNP \
	--variant $INDEL \
	-o 01-merged/$OUT.vcf \
	-genotypeMergeOptions Unsorted \
	-nt 24
 
 
done


##### merge all files together in one VCF

#java -Xmx10g -jar /data/ngs/bin/GATK_3.2-2/GenomeAnalysisTK.jar \
#-T CombineVariants \
#-K /data/ngs/bin/GATK_3.2-2/yong.li_uniklinik-freiburg.de.key -et NO_ET \
#-R /data/ngs/programs/bundle_2.8/ucsc.hg19.fasta \
#--variant 01-merged/20_9_1.vcf \
#--variant 01-merged/20_9_2.vcf \
#--variant 01-merged/20_9_3.vcf \
#--variant 01-merged/30_13_1.vcf \
#--variant 01-merged/30_13_2.vcf \
#--variant 01-merged/30_13_3.vcf \
#--variant 01-merged/498_2.vcf \
#--variant 01-merged/498_3.vcf \
#--variant 01-merged/498_4.vcf \
#--variant 01-merged/75_1.vcf \
#--variant 01-merged/75_2.vcf \
#--variant 01-merged/75_3.vcf \
#--variant 01-merged/86_24_1.vcf \
#--variant 01-merged/97_1_1.vcf \
#--variant 01-merged/97_1_2.vcf \
#--variant 01-merged/97_1_3.vcf \
#-o 2014-08-06-all.vcf \
#-genotypeMergeOptions Unsorted -nt 12
