
## recive in and out directories
inDir=$1
outDir=$2

mkdir -p $outDir



#define path
gemini="/data/ngs/bin/gemini/anaconda/bin/gemini"

#extract number of possible cpu
CPU=$(cat /proc/cpuinfo | grep processor | wc -l)


#for i in $(ls $1/*-vep-snpeff.vcf.gz)
#do

#        filename=$(basename "$i")
#        filename="${filename%.*}"
#        filename="${filename%.*}"
#       

#	echo $filename
#	echo $i

#	gemini load \
#	-t VEP \
#	--cores $CPU \
#	-v $i \
#	$outDir/$filename.db
#done



# old version didn't annotate cdna position due to snpEff 
for i in $(ls $1/*-vep-snpeff.vcf.gz)
do

        filename=$(basename "$i")
#        filename="${filename%.*}"
        filename="${filename%-vep-snpeff.vcf.gz}"
       
	echo $filename
	echo $i

	gemini load \
	-t VEP \
	--cores $CPU \
	-v $i \
	$outDir/$filename.db
done
