########  scheduler

	#wait a second for next job to start to avoid conflicts
	sleep 1

	#job controll, no more than maxJobs in paralell
	#check if jobs running is equal to max. if so sleep 5 sec and check again
	while [ $(jobs | wc -l) -ge $maxJobs ]
	do
		sleep 1
	done




######## read in file line wise and split to array

cat $bedFile | head | while IFS=$'\t' read -a chr
do


  echo ${chr[2]}

done
