mkdir -p 05-databases

#define path
gemini="/data/ngs/bin/gemini/bin/gemini"

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)


for i in $(ls $1/*-vep.vcf)
do

        filename=$(basename "$i")
        filename="${filename%.*}"

	echo $filename

	#add command
	#if ped-file exists
	if [ -e "ped/$filename.ped" ]
	then
		$gemini load \
		-t VEP \
		--cores $CPU \
		-v $i \
		-p ped/$filename.ped \
		05-databases/$filename.db
	else
		$gemini load \
		-t VEP \
		--cores $CPU \
		-v $i \
		05-databases/$filename.db
	fi
done
