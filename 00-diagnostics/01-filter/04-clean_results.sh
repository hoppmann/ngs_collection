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
		rm $i
	fi
done

# #### show line number of files
# ######## files with results"
# for i in $(ls $inputDir/*.out)
# do
#
# 	## write file name
# 	echo "######### $i"
#
# 	## write out number of lines of file
# 	numberlines=$( wc -l $i | awk '{ print $1 }' )
# 	echo "$numberlines lines in file (one line corresponds to header only)"
#
# 	## print out header
# 	echo -e "clinvar\tHGMD\tlineNumber\tclinvarPhenotype\tHGMDPhenotype"
#
# 	cat $i | awk -F $'\t' ' $33 ~ /(likely-)pathogenic/ || $41 ~ /DM/ {print "#### "$33 "\t" $41 " ####\t\t" $1 "\t" $36 "\t" $40}'
#
# 	echo ""
#
# done > $inputDir/summary.txt






#
# #### create single excel sheet of all _mod files
# sheets=""
# for i in $(ls $inputDir/*_mod.*)
# do
#
#         sheets="$sheets $i"
#
# done
# bn=$(basename $inputDir )
#
# if [ -e $inputDir/$bn.xlsx ]
# then
# 	rm $inputDir/$bn.xlsx
# fi
#
# ~/.local/bin/csv2xls \
# $sheets \
# -o $inputDir/$bn.xlsx \
# -d '	'





# awk -F $'\t' '{ print $33  "\t" $41 }'
