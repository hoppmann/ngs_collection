diag="/data/programs/scripts/hoppmann/00-diagnostics/"
panels="/data/programs/scripts/hoppmann/00-diagnostics/inSilicoPanels/"



for i in $(sqlite3 08-database/consensus-vep-snpEff.db "select * from samples" | cut -f 3 -d "|")
do
	if [[ $i == *"CVID"* ]]
	then
		$diag/00-filter.sh $i $panels/CCI/CVID/ exome &
	fi

	if [[ $i == *"HIE"* ]]
	then
		$diag/00-filter.sh $i $panels/CCI/HIES/ exome &
	fi

	if [[ $i == *"IBD"* ]]
	then
		$diag/00-filter.sh $i $panels/CCI/IBD/ exome &
	fi

	if [[ $i == *"CMC"* ]]
	then
		$diag/00-filter.sh $i $panels/CCI/CMC/ exome &
	fi

	if [[ $i == *"PID"* ]]
	then
		$diag/00-filter.sh $i $panels/CCI/PID/ exome &
	fi

done


wait

echo "Done"
