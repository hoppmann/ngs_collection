dir="01-vcf"


#mkdir -p $dir/oldReadgroup

#mv $dir/* $dir/oldReadgroup

for i in $(ls $dir/oldReadgroup/*freebayes*)
do
	file=$(basename $i)
	outFile="$dir/$file"
echo $outFile

	cat $i | awk '{ if ($1 ~ /^#CHROM/) { print ( $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" "freebayes" ) } else print } ' > $outFile
done


for i in $(ls $dir/oldReadgroup/*gatk*)
do

	file=$(basename $i)
	outFile="$dir/$file"
echo $outFile
	cat $i | awk '{ if ($1 ~ /^#CHROM/) { print ( $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" "gatk" ) } else print } ' > $outFile
done

for i in $(ls $dir/oldReadgroup/*platypus*)
do
	
	file=$(basename $i)
	outFile="$dir/$file"
echo $outFile
	cat $i | awk '{ if ($1 ~ /^#CHROM/) { print ( $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" "platypus" ) } else print } ' > $outFile
done

for i in $(ls $dir/oldReadgroup/*samtools*)
do
	
	file=$(basename $i)
	outFile="$dir/$file"	
	echo $outFile
	cat $i | awk '{ if ($1 ~ /^#CHROM/) { print ( $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" "samtools" ) } else print } ' > $outFile
done

