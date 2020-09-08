###################
########  scheduler
## define the maximum number of jobs to be executed in parallel
maxJobs=5
for i in {1..10}
do
	#wait a second for next job to start to avoid conflicts
	sleep 1

	## here add actual process to be executed and let it run in background "&"
	echo $i && sleep 20 &


	#job controll, no more than maxJobs in paralell
	#check if jobs running is equal to max. if so sleep 5 sec and check again
	while [ $(jobs | wc -l) -ge $maxJobs ]
	do
		sleep 5
	done

done



##################################################
######## read in file line wise and split to array

cat $fileIn | head | while IFS=$'\t' read -a chr
do

	echo ${chr[@]}
	echo ${chr[2]}

done



###########################################
######## get the path of the current script

SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPTPATH=`dirname $SCRIPT`
echo $SCRIPTPATH




#################################################
######## extract from match till frist emtpy line

awk '/Bad covered exons/,/^$/' SGA/ACAN.csv
sed -n '/Bad covered/,/^$/p' SGA/ACAN.csv
