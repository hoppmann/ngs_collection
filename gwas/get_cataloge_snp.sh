#!/bin/bash

if [ -z $1 ]
	then EXT=""
	else EXT="-$1"
fi
	
if [ -z $2 ]
	then EXT2=""
	else EXT2="-$2"
fi

FILES="adjusted unadjusted"


for FILE in $FILES
do
	echo "Processing..."
	echo gwama-ALL-${FILE}.out
	echo catalog_snps${EXT}${EXT2}-${FILE}.csv
	head -n 1 gwama-ALL-${FILE}.out > catalog_snps${EXT}${EXT2}-${FILE}.csv
	for i in $(cat snps)
	do
		grep -w $i gwama-ALL-${FILE}.out
	done >> catalog_snps${EXT}${EXT2}-${FILE}.csv
done
