mkdir -p 03-family

echo "family 14_1"
java -Xmx10g -jar /data/ngs/bin/GATK/3.2-2/GenomeAnalysisTK.jar \
-T CombineVariants -K \
/data/ngs/bin/GATK/3.2-2/yong.li_uniklinik-freiburg.de.key \
-et NO_ET \
-R /data/ngs/programs/bundle_2.8/ucsc.hg19.fasta \
--variant 00-vcf/21_261_1.vcf --variant 03-vcf/14_1.vcf --variant 03-vcf/4064.vcf \
-o 03-family/14_1.vcf \
-genotypeMergeOptions Unsorted \
-nt 24

