for i in $(squeue -u hoppmann -h -t PD -o %i)
do
	scontrol update jobid=$i partition="biom" NodeList="imbi4"
done
