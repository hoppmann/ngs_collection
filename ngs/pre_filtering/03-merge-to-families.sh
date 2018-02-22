mkdir -p 03-family

echo "family 719"
/dsk/localall/lib/java/jre1.8.0_131_x64/bin/java -Xmx10g -jar /data/ngs/bin/GATK/3.6.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /data/ngs/programs/bundle_2.8/ucsc.hg19.fasta \
--variant 01-torrentSuite_vcf/reRun/719_1.vcf \
--variant 01-torrentSuite_vcf/reRun/719_2.vcf \
--variant 01-torrentSuite_vcf/reRun/719_3.vcf \
-o 03-family/719.vcf \
-genotypeMergeOptions Unsorted \
-nt 20

