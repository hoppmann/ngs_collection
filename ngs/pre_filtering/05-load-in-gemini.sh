outDir=02-databases
mkdir -p $outDir

#define path
gemini="/data/ngs/bin/gemini/anaconda/bin/gemini"

##extract number of possible cpu
#CPU=$(cat /proc/cpuinfo | grep processor | wc -l)


#for i in $(ls $1/*-vep.vcf.gz	)
#do

#        filename=$(basename "$i")
#        filename="${filename%.*}"

#	echo $filename

#	#add command
#	#if ped-file exists
#	if [ -e "ped/$filename.ped" ]
#	then
#	echo "bitch"
##		$gemini load \
##		-t VEP \
##		--cores $CPU \
##		-v $i \
##		-p ped/$filename.ped \
##		05-databases/$filename.db
#	else
#		echo $i
#		$gemini load \
#		-t VEP \
#		--cores $CPU \
#		-v $i \
#		05-databases/$filename.db
#	fi
#done



#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)


for i in $(ls $1/*-vep-snpeff.vcf.gz	)
do

        filename=$(basename "$i")
        filename="${filename%.*}"
        filename="${filename%.*}"
       

	echo $filename
	echo $i

	gemini load \
	-t all \
	--cores $CPU \
	-v $i \
	$outDir/$filename.db
done
