#!/bin/bash

folderIn=$1
fileOut=$2



if [ $(ls $folderIn/* | wc -l ) -gt 1 ]
then
	rm  $fileOut
fi


for i in $(ls $folderIn/*.txt)
do
	sed -i 's/ //g' $i
	sed -i '/^\s*$/d' $i
done




if  [ $(ls $folderIn/* | wc -l ) -gt 1 ]
then
	for i in $(ls $folderIn/*.txt)
	do
		cat $i
	done > $folderIn/tmp

	cat $folderIn/tmp | sort | uniq > $fileOut
	rm $folderIn/tmp
fi
