#!/bin/bash


fileIn="$1"

dir=$(dirname "${fileIn}")
tmpFile=${dir}/tmp


omimDB="/data/public_resources/OMIM/2018-04/01-createDB/easyDB.db"


function join_by { local IFS="$1"; shift; echo "$*"; }


for curGene in $(cat $fileIn | cut -f 1)
do

	allmoi=$(sqlite3 $omimDB "select AR, AD, XLR, XLD from omim where gene == '$curGene'")
	IFS="|" read -r -a moiSplit <<< "$allmoi"


	moiColect=()
	if [ ! -z ${moiSplit[0]} ]
	then
		moiColect+=("AR")
	fi

	if [ ! -z ${moiSplit[1]} ]
	then
		moiColect+=("AD")
	fi

	if [ ! -z ${moiSplit[2]} ]
	then
		moiColect+=("XLR")
	fi

	if [ ! -z ${moiSplit[3]} ]
	then
		moiColect+=("XLD")
	fi

	moi=$(join_by "," ${moiColect[@]})
	echo -e "$curGene\t$moi"

done > $tmpFile

mv $tmpFile $fileIn
chmod 755 $fileIn
