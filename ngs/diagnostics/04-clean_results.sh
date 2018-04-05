clear
inputDir="$1"

for i in $(ls $inputDir/*)
do

	length=$(wc -l $i | awk '{ print $1 }')
	if [ $length -eq 1 ]
	then
		echo $i
		# rm $i
	fi
done


#### show line number of files
echo ""
echo "######## files with results"
wc -l $inputDir/*
echo ""
echo ""
