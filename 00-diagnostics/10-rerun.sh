clear


for i in {38..42}
do
	folder=$(ls -d *Run$i/)
	cd $folder
	# ls 00-fastq
	# cd ..


	dir=$(ls -d NextSeq_Run[0-9][0-9]/)
	# dir=$(ls -d NextSeq_Run[1-9]/)
	# echo $dir
	# ls $dir
	bn=$(basename $dir)
	time=$((${bn##*Run}-18))
	echo "$fn"
	echo -e "copying alamut to $dir/07-annotation/"

	# echo $time

	mkdir -p  $dir/07-annotation

	cp -r /tempdata/ge/temp/hoppmann/NextGP/$dir/08-annotation/Alamut/ \
	$dir/07-annotation/

	echo -e "preparing new pipeline\n"
	rm -r pipeline/
	rm -r slurmLog/
	cp /data/programs/scripts/hoppmann/00-diagnostics/01-run_analysis.sh .
	. 01-run_analysis.sh pipe

	sbatch --begin=now+${time}hour pipeline/01-master.sh

	cd ..

done
