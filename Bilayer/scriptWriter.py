import os
import numpy as np
import pdb
import random
import sys

"""
scriptWriter class 
stores information about the filename
has methods to write submission scipts (initial, continue, and repeat files) 
for Cori and Rahman clusters
"""

class scriptWriter():
    def __init__(self, filename):
        self._filename = filename

    # Cori writing
    def write_Cori_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STsbatch.sbatch'),'w')
        elif MDrun:
            init_file = open((filename + 'MDsbatch.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write('#!/bin/bash -l\n')
        init_file.write('#SBATCH -p low\n')
        init_file.write('#SBATCH -N 8\n')
        if STrun:
            init_file.write('#SBATCH -t 07:00:00\n')
        elif MDrun:
            init_file.write('#SBATCH -t 09:00:00\n')
        else:
            pass
        init_file.write('#SBATCH -J {}\n'.format(filename))
        init_file.write('#SBATCH -o {}\n'.format(filename))
        init_file.write('#SBATCH -C haswell\n')
        init_file.write('#SBATCH --mail-type=ALL\n')
        init_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        init_file.write('module load gromacs/5.1.2\n')
        if STrun:
            init_file.write('srun -n 128 -c 4 mdrun_mpi_sp -deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
        elif MDrun:
            init_file.write('srun -n 128 -c 4 mdrun_mpi_sp -deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
        else:
            pass
            
        init_file.close()
        
        # Write continue script
        if STrun:
            cont_file = open((filename + 'STcont.sbatch'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDcont.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write('#!/bin/bash -l\n')
        cont_file.write('#SBATCH -p low\n')
        cont_file.write('#SBATCH -N 8\n')
        if STrun:
            cont_file.write('#SBATCH -t 07:00:00\n')
        elif MDrun:
            cont_file.write('#SBATCH -t 09:00:00\n')
        else:
            pass
        cont_file.write('#SBATCH -J {}\n'.format(filename))
        cont_file.write('#SBATCH -o {}\n'.format(filename))
        cont_file.write('#SBATCH --mail-type=ALL\n')
        cont_file.write('#SBATCH -C haswell\n')
        cont_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        cont_file.write('module load gromacs/5.1.2\n')
        if STrun:
            cont_file.write('srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/{}/ST_{}.cpt \\\n'.format(filename,filename))
            cont_file.write('-s $SCRATCH/Trajectories/{}/ST_{}.tpr \\\n'.format(filename,filename))
            cont_file.write('-deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
        elif MDrun:
            cont_file.write('srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/{}/md_{}.cpt \\\n'.format(filename,filename))
            cont_file.write('-s $SCRATCH/Trajectories/{}/md_{}.tpr \\\n'.format(filename,filename))
            cont_file.write('-deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
        else:
            pass
        
        cont_file.close()
        
        
        # Write repeat script
        if STrun:
            repeat_file = open((filename + 'STrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}STsbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}STcont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}STcont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        elif MDrun:
            repeat_file = open((filename + 'MDrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}MDsbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}MDcont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}MDcont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        else:
            pass
        repeat_file.close()

    def write_Rahman_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STpbs.pbs'),'w')
        elif MDrun:
            init_file = open((filename + 'MDpbs.pbs'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write("#!/bin/sh -l\n")
        init_file.write("#PBS -N {}\n".format(filename))
        init_file.write("#PBS -l nodes=1:ppn=16\n")
        if STrun:
            init_file.write("#PBS -l walltime=14:00:00\n")
        elif MDrun:
            init_file.write("#PBS -l walltime=24:00:00\n")
        else: 
            pass
        init_file.write("#PBS -q low\n")
        init_file.write("#PBS -m abe\n")
        init_file.write("#PBS -M alexander.h.yang@vanderbilt.edu\n")
        init_file.write("cd $PBS_O_WORKDIR\n")
        init_file.write("echo `cat $PBS_NODEFILE`\n")
        init_file.write("module load gromacs/5.1.0\n")
        init_file.write("cd ~/Trajectories/{}\n".format(filename))
        if STrun:
            init_file.write("gmx mdrun -ntomp 8 -gpu_id 0 -deffnm ST_{} >& out.log\n".format(filename))
        elif MDrun:
            init_file.write("gmx mdrun -ntomp 8 -gpu_id 0 -deffnm md_{} >& out.log\n".format(filename))
        else:
            pass
        init_file.close()

        if STrun:
            cont_file = open((filename + 'STcont.pbs'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDcont.pbs'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write("#!/bin/sh -l \n")
        cont_file.write("#PBS -N {}\n".format(filename))
        cont_file.write("#PBS -l nodes=1:ppn=16\n")
        cont_file.write("#PBS -l walltime=04:00:00\n")
        cont_file.write("#PBS -q low\n")
        cont_file.write("cd $PBS_O_WORKDIR\n")
        cont_file.write("echo `cat $PBS_NODEFILE`\n")
        cont_file.write("module load gromacs/5.1.0\n")
        cont_file.write("cd ~/Trajectories/{}\n".format(filename))
        cont_file.write("gmx mdrun -ntomp 8 -gpu_id 0 -append \ \n")
        if STrun:
            cont_file.write("-s ST_{}.tpr \ \n".format(filename))
            cont_file.write("-cpi ST_{}.cpt \ \n".format(filename))
            cont_file.write("-deffnm ST_{}\n".format(filename))
        elif MDrun:
            cont_file.write("-s md_{}.tpr \ \n".format(filename))
            cont_file.write("-cpi md_{}.cpt \ \n".format(filename))
            cont_file.write("-deffnm md_{}\n".format(filename))
        else:
            pass
        cont_file.close()
