export item=`sbatch JOBNAMEsbatch.sbatch`
#export item=`sbatch --dependency=afterany:2418639 JOBNAMEcont.sbatch`
for i in {0..2}
do
	item=$(sbatch --dependency=afterany:${item:20:7} JOBNAMEcont.sbatch)
done
