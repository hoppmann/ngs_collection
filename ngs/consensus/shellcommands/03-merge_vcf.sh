dirIn="01-vcf"
dirOut="02-merged"

mkdir -p $dirOut


for i in $(ls $dirIn/*.vcf)
do

	file=$(basename $i)
	curPat=${file%_*.vcf}
	
	srun -p genepi \
	/dsk/localall/lib/java/jre1.8.0_111_x64/bin/java \
	-jar /data/ngs/bin/GATK/3.6.0/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
	--variant $dirIn/$curPat\_freebayes.vcf \
	--variant $dirIn/$curPat\_gatk.vcf \
	--variant $dirIn/$curPat\_platypus.vcf \
	--variant $dirIn/$curPat\_samtools.vcf \
	-genotypeMergeOptions Unsorted \
	-o $dirOut/$curPat.vcf &



#/data/ngs/bin/vcftools/0.1.12b/bin/vcf-merge $dirIn/FG110_freebayes.vcf  $dirIn/FG110_gatk.vcf  $dirIn/FG110_platypus.vcf  $dirIn/FG110_samtools.vcf | bgzip -c > out.vcf.gz



done
