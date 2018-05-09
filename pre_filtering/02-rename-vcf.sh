#replace Unknown with actual patient ID for all vcf

pathToVCF=$1

#create out dir if not allready exist
mkdir -p 02-renamed

for i in $(ls $pathToVCF/*.vcf)
do 
	echo $i
	filename=$(basename "$i")
	filename="${filename%.*}"
	echo $filename
	cp $i 02-renamed/$filename.vcf
	sed -i "s/Unknown/${filename}/g" 02-renamed/$filename.vcf
done
