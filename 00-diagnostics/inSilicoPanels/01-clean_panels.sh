#!/bin/bash

folderIn=$1
fileOut=$2



if [ $(ls $folderIn/* | wc -l ) -gt 1 ]
then
	rm  $fileOut
fi

echo ""
echo "######## $folderIn ########"

for i in $(ls $folderIn/*.txt)
do
	echo ""
	echo "Cleaning $i"

	sed -i 's/ //g' $i
	sed -i '/^\s*$/d' $i
	sed -i 's/?//g' $i
	sed -i 's/+/,/g' $i

	cat $i | sort | uniq > tmp
	mv tmp $i

	echo "Annotating genes"
	. 02-annotate_panels.sh $i

done




if  [ $(ls $folderIn/* | wc -l ) -gt 1 ]
then
	for i in $(ls $folderIn/*.txt)
	do
		cat $i
	done > $folderIn/tmp

	cat $folderIn/tmp | sort | uniq > $fileOut
	rm $folderIn/tmp

	for i in $(ls $folderIn/*.txt)
	do
		dos2unix $i
	done


fi
