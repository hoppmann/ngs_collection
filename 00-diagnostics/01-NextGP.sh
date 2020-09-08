#!/bin/bash
#SBATCH -c 24
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6



######## general variables
maxCPU=6
totalCPU=20
# inDir="NextSeq_Run1_171129_NB500913_0003_AH3N3TBGX5"
inDir=$(ls -d *_Run*_[0-1]*/)
outDir="00-fastq"
NGPFolder=${inDir%%_[0-9]*}
#### if folder does not exist take generic result dir
if [ -z $NGPFolder ]
then
	NGPFolder=$(basename $(pwd))
fi

# bcl2fastq="/data/programs/bin/ngs/bcl2fastq/2.20.0/bin/bcl2fastq"
bcl2fastq="/data/programs/bin/ngs/bcl2fastq/2.20.0.422-Deb10-vanilla/bin/bcl2fastq"

bedFileV4="/data/public_resources/bed/SureSelect_V4_S03723314_Regions.bed"
bedFileV6="/data/public_resources/bed/SureSelectExomeV6_S07604514_Padded.bed"
bedFileV7="/data/public_resources/bed/SureSelectExomeV7_S31285117_Padded_hg19.bed"
bedFileTwist="/data/public_resources/bed/Twist_Exome_RefSeq_targets_hg19_extended_by_30bp.bed"

# customBedFile="/data/studies/01_ngs/01_Kinderklinik/Miriam/00-Iran_cohort/00-bed/S07604514_v6_bed/S07604514_Padded.bed"


# customBedFile="/data/public_resources/bed/SureSelectExomeV7_S31285117_Padded_hg19.bed"
customBedFile=$bedFileV6
# customBedFile=/data/studies/01_ngs/01_Kinderklinik/misc/2020-05-28_BGI/hg19_Agilent_V6_region.sort.bed
#### CCI BED Files
# customBedFile="../00-bed/hg19/CCI_PID_v5_3195801_Covered.hg19.bed"
# customBedFile="/data/studies/01_ngs/02_CCI/00-bed/hg19/CCI_PID_v5_3195801_Covered.hg19.bed"

#### BGI BED Files
# customBedFile="/data/studies/01_ngs/01_Kinderklinik/misc/2020-05-28_BGI/hg19_Agilent_V6_region.sort.bed"


if [[ $NGPFolder == *"V7"* ]]
then
	bedFile=$bedFileV7
elif [[ $NGPFolder == *"V6"* ]]
then
	bedFile=$bedFileV6
elif [[ $NGPFolder == *"V4"* ]]
then
	bedFile=$bedFileV4
elif [[ $NGPFolder == *"TWIST"* ]]
then
	bedFile=$bedFileTwist
else
	bedFile=$customBedFile
fi


##############################
####### run bcl2fastq ########
##############################



if [[ "$*" == *"fastq"* ]]
then

	#### retrieve sequence size for bcl2fastq index extraction
	read1=$(cat "$inDir/RunInfo.xml" | grep "Read Number=\"1\"")
	read1=${read1##*NumCycles=\"}
	read1=${read1%%\"*}

	index1=$(cat "$inDir/RunInfo.xml" | grep "Read Number=\"2\"")
	index1=${index1##*NumCycles=\"}
	index1=${index1%%\"*}


	index2=$(cat "$inDir/RunInfo.xml" | grep "Read Number=\"3\"")
	index2=${index2##*NumCycles=\"}
	index2=${index2%%\"*}

	read2=$(cat "$inDir/RunInfo.xml" | grep "Read Number=\"4\"")
	read2=${read2##*NumCycles=\"}
	read2=${read2%%\"*}


	# echo "$read1 $index1 $index2 $read2"

	#### run demultiplexing

	if [[ $NGPFolder == *"TWIST"* ]]
	then

		echo "Running TWIST bcl2fastq"

		$bcl2fastq \
		-R $inDir \
		-o $outDir \
		-r $maxCPU \
		-p $maxCPU \
		-w $maxCPU \
		--create-fastq-for-index-reads \
		--no-lane-splitting \
		--mask-short-adapter-reads 8 \
		--use-bases-mask y${read1},i${index1},i${index2},y${read2} \
		--ignore-missing-bcls

	else

		echo "Running basic bcl2fastq"

		$bcl2fastq \
		-R $inDir \
		-o $outDir \
		-r $maxCPU \
		-p $maxCPU \
		-w $maxCPU \
		--create-fastq-for-index-reads \
		--no-lane-splitting \
		--mask-short-adapter-reads 8 \
		--use-bases-mask y${read1},i${index1},y${index2},y${read2} \
		--ignore-missing-bcls


		#### rename files
		# I1 -> I1
		# R1 -> R1
		# R2 -> I2
		# R3 -> R2

		for i in $( ls $outDir/*R2*.fastq.gz)
		do
			R2=$i
			I2=${R2/R2/I2}
			R3=${R2/R2/R3}

			mv $R2 $I2
			mv $R3 $R2

		done

	fi
fi









 ################################
 ######## execute NextGP ########
 ################################

 if [[ "$@" == *"pipe"* ]]
 then

 ######## prepare fastq.list

 for i in $(ls 00-fastq/*R1*)
 do
   R1=$i
   R2=${R1/R1/R2}
   bn=$(basename $R1)
   fn=${bn%*_S*}
   fn=${fn%*_R*}


   if [[ $fn != "Undetermined" ]]
   then
     echo -e "$fn\t$R1\t$R2"
   fi

done > fastq.list





#############################################
####### run NextGP from fastq or BAM ########
#############################################

 #		first		Number of first pipeline step to start from. Default: 01.
 #				Possible steps:
 #					01: Preprocess
 #					02: alignment;
 #					03: remove duplicates;
 #					04: indel realignment;
 #					05: base recalibration;
 #					06: metrices;
 #					07: variant calling;
 #					08: annotation;
 #					09: database generation
 #					10: annotate Alamut and Post Gemini annotations
 #					11: filter variants
 #		last		Last step of Pipeline to run. Default 10.
 #		noMail		If set, no notice of failing or finishing will be send to the corresponding e-mail adress deposed in slurm.
 #		skip		Skips running alamut batch. Only use if an old annotation file is available.
 #		after		":"-seperated list of Slurm PIDs for afterok-dependency for the 05-Job.
 #

 	# # ####### bamList setting ########
 	# /data/programs/pipelines/NextGP/NextGP \
 	# -bedFile $bedFile \
 	# -outDir $NGPFolder \
 	# -first 1 \
 	# -last 10 \
 	# -bamList bam.list \
 	# -exome

 	####### NextSeq setting ########
 	/data/programs/pipelines/NextGP/NextGP \
 	-bedFile $bedFile \
 	-outDir $NGPFolder \
 	-first 1 \
 	-last 11 \
 	-fastqList fastq.list \
 	-exome \
	# -skip
	# -slurmPartition biom \
 	# -restrict imbi3 \
	# -skip \
 	# -restrict imbi6 \
 	# -exclude imbi8 \


 fi




############################################
######## execute NextGP from VCF on ########
############################################

if [[ "$@" == *"vcf"* ]]
then

	######## prepare fastq.list

	for curVcf in $(ls 00-vcf/*.vcf)
	do
	  bn=$(basename $curVcf)
	  fn=${bn%*.vcf}

	  echo "$fn	$curVcf"

  done > vcf.list

	/data/programs/pipelines/NextGP/NextGP \
	-outDir $NGPFolder \
	-first 08 \
	-last 11 \
	-exome \
	-vcfList vcf.list \
	# -slurmPartition biom \
 	# -restrict imbi12 \
	# -skip \
	# -exclude imbi5 \
	# -exclude imbi11 \



fi




























 ################################
 ######## start pipeline ########
 ################################

 if [[ "$@" == *"run"* ]]
 then
 	. pipeline/01-master.sh
 fi

























echo "Script finished!"
