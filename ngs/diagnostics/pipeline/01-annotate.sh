#!/bin/bash
#SBATCH --mailType=FAIL
#SBATCH --cpus-per-task=6


## recieve in and out dir
inDir=$1
OUT=$2


# progs
#vep="/data/ngs/bin/VEP/77/variant_effect_predictor.pl"
vep="/data/ngs/bin/VEP/89/vep"
cacheDir="/data/ngs/bin/VEP/cache"
vt="/data/ngs/bin/vt-master/vt"

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)
CPU=6

mkdir -p $OUT

for i in $(ls $inDir/*.vcf)
do




    	echo $i
	filename=$(basename "$i")
        filename="${filename%.*}"

	# tabix files
	bgzip -c $i > $i.gz
	tabix $i.gz



        #prepare for anotation-c
	echo $i
       zless $i \
       | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
       | $vt decompose -s - \
       | $vt normalize -r /data/ngs/resources/bundle_2.8/ucsc.hg19.fasta - > $OUT/$filename.vt.vcf


        # annotate using VEP
	$vep \
	-i $OUT/$filename.vt.vcf \
	-o $OUT/$filename-vep.vcf \
	--fork $CPU \
	--cache \
	--merged \
	--offline \
	--dir_cache $cacheDir \
	--cache_version 89 \
	--assembly GRCh37 \
	--plugin dbNSFP,/data/ngs/resources/dbNSFP/2.9.3/dbNSFP2.9.3_hg19.gz,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,MetaSVM_pred,MetaLR_pred,PROVEAN_pred,M-CAP_pred,REVEL_score,clinvar_clnsig,clinvar_trait \
	--force_overwrite \
	--sift b \
	--polyphen b \
	--symbol \
	--numbers \
	--biotype \
	--total_length \
	--check_existing \
	--canonical \
	--pubmed \
	--vcf



        ## annotate using snpEff
	echo "Running SNPEff!"
	java -jar /data/ngs/bin/snpEff/4.3p/snpEff.jar \
	GRCh37.75 \
	-classic \
	-formatEff \
	-noStats \
	$OUT/$filename-vep.vcf \
	> $OUT/$filename-vep-snpeff.vcf

	# indexing and ziping out file
	echo "Indexing files!"
	bgzip -c $OUT/$filename-vep-snpeff.vcf > $OUT/$filename-vep-snpeff.vcf.gz
	tabix -p vcf $OUT/$filename-vep-snpeff.vcf.gz








done
