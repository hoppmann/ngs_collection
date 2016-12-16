clear

dir="03-consensus"

for i in $(ls $dir/oldReadgroups/*)
do

	file=$(basename $i)
	curPat=${file%.*}
	outFile="$dir/$file"
	echo "Processing $file"


	cat $i | awk -v curPat="$curPat" '{ if ($1 ~ /^#CHROM/) { print ( $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" curPat ) } else print } ' > $outFile


done
