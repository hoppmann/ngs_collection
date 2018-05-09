#clear
inputDir="$1"

if [ -e $inputDir/summary.txt ]
then
	rm $inputDir/summary.txt
fi

for i in $(ls $inputDir/*)
do

	length=$(wc -l $i | awk '{ print $1 }')
	if [ $length -eq 1 ]
	then
		# echo $i
		rm $i
	fi
done

#### show line number of files
# echo ""
# echo "######## files with results"
for i in $(ls $inputDir/*)
do

	## write file name
	echo "######### $i"

	## write out number of lines of file
	numberlines=$( wc -l $i | awk '{ print $1 }' )
	echo "$numberlines lines in file (one line corresponds to header only)"

	## print out header
	echo -e "clinvar\tHGMD\tlineNumber\tclinvarPhenotype\tHGMDPhenotype"

	cat $i | awk -F $'\t' ' $33 ~ /(likely-)pathogenic/ || $41 ~ /DM/ {print "#### "$33 "\t" $41 " ####\t\t" $1 "\t" $36 "\t" $40}'
	# cat $i | awk -F $'\t' ' $33 ~ /pathogenic/ {print $33 "\t" $41 }'

	echo ""

done > $inputDir/summary.txt


# awk -F $'\t' '{ print $33  "\t" $41 }'
