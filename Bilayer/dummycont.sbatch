#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 8
#SBATCH -t 10:00:00
#SBATCH -J JOBNAME
#SBATCH -o JOBNAME
#SBATCH --mail-type=ALL
#SBATCH -C haswell
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
module load gromacs/5.1.2
srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/JOBNAME/md_JOBNAME.cpt \
-s $SCRATCH/Trajectories/JOBNAME/md_JOBNAME.tpr \
-deffnm $SCRATCH/Trajectories/JOBNAME/md_JOBNAME >& out.log
