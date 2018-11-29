#!/bin/bash
#SBATCH -c 24
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6



###################################
######## general variables ########
###################################

maxCPU=6
totalCPU=24
inDir="NextSeq_Run14_180522_NB500913_0016_AHWNHKBGX5"
NGPFolder=${inDir%%_[0-9]*}
outDir="00-fastq"




###############################
######## run bcl2fastq ########
###############################

/data/programs/bin/ngs/bcl2fastq/2.20.0/bin/bcl2fastq \
-R $inDir \
-o $outDir \
-r $maxCPU \
-p $maxCPU \
-w $maxCPU \
--create-fastq-for-index-reads \
--no-lane-splitting



####################################
######## prepare fastq.list ########
####################################

for i in $(ls 00-fastq/*R1*)
do
  R1=$i
  R2=${R1/R1/R2}
  bn=$(basename $R1)
  fn=${bn%*_S*}


  if [[ $fn != "Undetermined" ]]
  then
    echo -e "$fn\t$R1\t$R2"
  fi

done > fastq.list

# ##################################
# ######## prepare bam.list ########
# ##################################
#
# for i in $(ls 00-bam/*.bam)
# do
#   BAM=$i
#   bn=$(basename $BAM)
#   fn=${bn%*.bam}
#
#
#   if [[ $fn != "Undetermined" ]]
#   then
#     echo -e "$fn\t$BAM"
#   fi
#
# done > bam.list

###########################
####### run NextGP ########
###########################

/data/programs/pipelines/NextGP/NextGP \
-bedFile /data/public_resources/bed/SureSelectExomeV6_S07604514_Padded.bed \
-outDir $NGPFolder \
-fastqList fastq.list \
-exclude fdm203 \
-solid \
-mail \
-exon

##############################q
######## run pipeline ########
##############################

/. pipeline/01-master.sh
