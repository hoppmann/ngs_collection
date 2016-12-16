folderIn=$1
folderIn="out/04-baseRecalibration"

outDir="01-vcf"
mkdir -p $outDir

log="log"
mkdir -p $log


#for i in $(ls $folderIn/FG110.bam)
for i in $(ls $folderIn/*.bam)	
do

logOut=$(basename $i)
logOut=${logOut%.*}
curPat=$logOut


	## call variants using samtools
	srun -p genepi \
	/data/ngs/bin/samtools/1.3.1/bin/samtools mpileup \
	--uncompressed \
	--VCF \
	--fasta-ref /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
	--positions /data/ngs/rohdaten/kinderklinik/misc/2016-08-17-Stoffwechsel_Panel/00-bed/04818-1453471441_Regions.bed \
	$i \
	| bcftools call -f GQ -vmO z -o $outDir/$logOut\_samtools.vcf.gz  > $log/$logOut.samtools.log & 
	
	

	## call variants using platypus and merge clustered variants
	srun -p genepi \
	python /data/ngs/bin/platypus/bin/Platypus.py \
	callVariants \
	--bamFiles=$i \
	--refFile=/data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
	--logFileName=platypus.log \
	--nCPU=2 \
	--mergeClusteredVariants=1 \
	--output=$outDir/$logOut\_platypus.vcf > $log/$logOut.platypus.log &


	## extract individuals from pipelinerun using gatk
	srun -p genepi \
	/data/ngs/bin/vcftools/0.1.12b/bin/vcftools \
	--recode \
	--indv $curPat \
	--vcf out/07-variantFiltering/combined_filtered.recode.vcf \
	--out $outDir/$logOut\_gatk.vcf > $log/$logOut.gatk.log &



	## freebayes
	srun -p genepi \
	/data/ngs/bin/freebayes/bin/freebayes \
	-f /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta \
	-C 5 \
	--genotype-qualities \
	$i \
	| /data/ngs/bin/freebayes/vcflib/bin/vcffilter \
	-f "QUAL > 20" \
	> $outDir/$logOut\_freebayes.vcf  &
	
	
done


### clean platypus from filtered reads

#for i in $(ls $outDir/*platypus*)
#do
#/data/ngs/bin/vcftools/0.1.12b/bin/vcftools \
#--vcf $i \
#--remove-filtered-all \
#--recode \
#--out $i.filtered
#done

## unzip all zipped files
#for i in $(ls $outDir/*.gz)
#do
#	bgzip -d $i
#done


## clean up stuff
#rename -f 's/.filtered.recode.vcf//' $outDir/*
#rename -f 's/.recode.vcf//' $outDir/*
#rm $outDir/*.log
#chmod 750 $outDir/*





