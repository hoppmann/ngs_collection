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

cat $bedFile | head | while IFS=$'\t' read -a chr
do

	echo ${chr{@}}
	echo ${chr[2]}

done
