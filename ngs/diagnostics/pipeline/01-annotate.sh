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
	--plugin dbNSFP,/data/ngs/resources/dbNSFP/2.9.3/dbNSFP2.9.3_hg19.gz,cds_strand,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,VEST3_score,VEST3_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore \
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
