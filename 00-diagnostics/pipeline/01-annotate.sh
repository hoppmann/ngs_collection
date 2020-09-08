#!/bin/bash
#SBATCH --mailType=FAIL
#SBATCH --cpus-per-task=6


## recieve in and out dir
inDir=$1
OUT=$2


# progs
vep="/data/programs/bin/ngs/VEP/96/vep"
cacheDir="/data/programs/bin/ngs/VEP/cache"
vt="/data/programs/bin/ngs/vt-master/vt"
snpEff="/data/programs/bin/ngs/snpEff/4.3t/snpEff.jar"
cacheVersion="96"
assembly="GRCh37"
species="homo_sapiens"


java8="/dsk/localall/lib/java/jre1.8.0_131_x64/bin/java"


#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)
CPU=12

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
	| $vt normalize -r /data/public_resources/bundle/2.8/ucsc.hg19.fasta - > $OUT/$filename.vt.vcf


        # annotate using VEP
	echo "#### Running VEP ####"
	echo ""
	$vep \
	-i $OUT/$filename.vt.vcf \
	-o $OUT/$filename-vep.vcf \
	--fork $CPU \
	--cache \
	--merged \
	--offline \
	--dir_cache $cacheDir \
	--cache_version $cacheVersion \
	--assembly $assembly \
	--species $species \
	--plugin dbNSFP,/data/public_resources/dbNSFP/4.0b1a/dbNSFP_hg19.gz,SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MetaLR_score,MetaLR_pred,MetaSVM_score,MetaSVM_pred,REVEL_score,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred \
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
#	--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE



    ## annotate using snpEff
	echo "#### Running SNPEff! ####"
	echo ""
	$java8 -jar $snpEff \
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

    #bgzip -c $OUT/$filename-vep.vcf > $OUT/$filename-vep.vcf.gz
    #tabix -p vcf $OUT/$filename-vep.vcf.gz






done


#cds_strand,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,VEST3_score,VEST3_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore
