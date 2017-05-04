import os
import numpy as np
import pdb
import random
import sys

"""
scriptWriter class 
stores information about the filename
has methods to write submission scipts (initial, continue, and repeat files) 
for Cori, Titan, and Rahman clusters
"""

class scriptWriter():
    def __init__(self, filename):
        self._filename = filename

    # Cori writing
    def write_Cori_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STCorisbatch.sbatch'),'w')
        elif MDrun:
            init_file = open((filename + 'MDCorisbatch.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write('#!/bin/bash -l\n')
        init_file.write('#SBATCH -p regular\n')
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
            cont_file = open((filename + 'STCoricont.sbatch'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDCoricont.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write('#!/bin/bash -l\n')
        cont_file.write('#SBATCH -p regular\n')
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
            repeat_file = open((filename + 'STCorirepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}STCorisbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}STCoricont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}STCoricont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        elif MDrun:
            repeat_file = open((filename + 'MDCorirepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}MDCorisbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}MDCoricont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}MDCoricont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        else:
            pass
        repeat_file.close()
    
    def write_Edison_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STEdisonsbatch.sbatch'),'w')
        elif MDrun:
            init_file = open((filename + 'MDEdisonsbatch.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write('#!/bin/bash -l\n')
        init_file.write('#SBATCH -p regular\n')
        init_file.write('#SBATCH -N 2\n')
        if STrun:
            init_file.write('#SBATCH -t 18:00:00\n')
        elif MDrun:
            init_file.write('#SBATCH -t 18:00:00\n')
        else:
            pass
        init_file.write('#SBATCH -J {}\n'.format(filename))
        init_file.write('#SBATCH --mail-type=ALL\n')
        init_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        init_file.write('#SBATCH -L SCRATCH\n')
        init_file.write('module load gromacs/5.1.2\n')
        init_file.write('export OMP_NUM_THREADS=1\n')
        if STrun:
            init_file.write('srun -n 96 -c 1 mdrun_mpi_sp -ntomp 1 -deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
        elif MDrun:
            init_file.write('srun -n 96 -c 1 mdrun_mpi_sp -ntomp 1 -deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
        else:
            pass
            
        init_file.close()
        
        # Write continue script
        if STrun:
            cont_file = open((filename + 'STEdisoncont.sbatch'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDEdisoncont.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write('#!/bin/bash -l\n')
        cont_file.write('#SBATCH -p regular\n')
        cont_file.write('#SBATCH -N 2\n')
        if STrun:
            cont_file.write('#SBATCH -t 07:00:00\n')
        elif MDrun:
            cont_file.write('#SBATCH -t 09:00:00\n')
        else:
            pass
        cont_file.write('#SBATCH -J {}\n'.format(filename))
        cont_file.write('#SBATCH --mail-type=ALL\n')
        cont_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        cont_file.write('#SBATCH -L SCRATCH\n')
        cont_file.write('module load gromacs/5.1.2\n')
        cont_file.write('export OMP_NUM_THREADS=1\n')
        if STrun:
            cont_file.write('srun -n 96 -c 1 mdrun_mpi_sp -ntomp 1 -append -cpi $SCRATCH/Trajectories/{}/ST_{}.cpt \\\n'.format(filename,filename))
            cont_file.write('-s $SCRATCH/Trajectories/{}/ST_{}.tpr \\\n'.format(filename,filename))
            cont_file.write('-deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
        elif MDrun:
            cont_file.write('srun -n 96 -c 1 mdrun_mpi_sp -ntomp 1 -append -cpi $SCRATCH/Trajectories/{}/md_{}.cpt \\\n'.format(filename,filename))
            cont_file.write('-s $SCRATCH/Trajectories/{}/md_{}.tpr \\\n'.format(filename,filename))
            cont_file.write('-deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
        else:
            pass
        
        cont_file.close()
        
        
        # Write repeat script
        if STrun:
            repeat_file = open((filename + 'STEdisonrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}STEdisonsbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}STEdisoncont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}STEdisoncont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        elif MDrun:
            repeat_file = open((filename + 'MDEdisonrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}MDEdisonsbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}MDEdisoncont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}MDEdisoncont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        else:
            pass
        repeat_file.close()


    def write_Rahman_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STRahmanpbs.pbs'),'w')
        elif MDrun:
            init_file = open((filename + 'MDRahmanpbs.pbs'),'w')
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
            cont_file = open((filename + 'STRahmancont.pbs'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDRahmancont.pbs'),'w')
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

    def write_Titan_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STTitanpbs.pbs'),'w')
        elif MDrun:
            init_file = open((filename + 'MDTitanpbs.pbs'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write("#!/bin/sh -l\n")
        init_file.write("#PBS -N {}\n".format(filename))
        init_file.write("#PBS -l nodes=8\n")
        if STrun:
            init_file.write("#PBS -l walltime=02:00:00\n")
        elif MDrun:
            init_file.write("#PBS -l walltime=02:00:00\n")
        else: 
            pass
        init_file.write("#PBS -j oe\n")
        init_file.write("#PBS -A MAT149\n")
        init_file.write("#PBS -m abe\n")
        init_file.write("#PBS -M alexander.h.yang@vanderbilt.edu\n")
        init_file.write("cd $PBS_O_WORKDIR\n")
        init_file.write("echo `cat $PBS_NODEFILE`\n")
        init_file.write("module load gromacs/5.1.0\n")
        init_file.write("export CRAY_CUDA_MPS=1\n")
        init_file.write("export OMP_NUM_THREADS=2\n")
        init_file.write("cd $MEMBERWORK/mat149/Trajectories/{}\n".format(filename))
        if STrun:
            init_file.write("aprun -n 64 -N 8 gmx_mpi mdrun -gpu_id 00000000 -deffnm ST_{} >& out.log\n".format(filename))
        elif MDrun:
            init_file.write("aprun -n 64 -N 8 gmx_mpi mdrun -gpu_id 00000000 -deffnm md_{} >& out.log\n".format(filename))
        else:
            pass
        init_file.close()

        if STrun:
            cont_file = open((filename + 'STTitancont.pbs'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDTitancont.pbs'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write("#!/bin/sh -l \n")
        cont_file.write("#PBS -N {}\n".format(filename))
        cont_file.write("#PBS -l nodes=8\n")
        cont_file.write("#PBS -l walltime=02:00:00\n")
        cont_file.write("#PBS -j oe\n")
        cont_file.write("#PBS -A MAT149\n")
        cont_file.write("#PBS -M alexander.h.yang@vanderbilt.edu\n")
        cont_file.write("cd $PBS_O_WORKDIR\n")
        cont_file.write("echo `cat $PBS_NODEFILE`\n")
        cont_file.write("module load gromacs/5.1.0\n")
        cont_file.write("export CRAY_CUDA_MPS=1\n")
        cont_file.write("export OMP_NUM_THREADS=2\n")
        cont_file.write("cd $MEMBERWORK/mat149/Trajectories/{}\n".format(filename))
        cont_file.write("aprun -N 8 -n 64 gmx_mpi mdrun -gpu_id 00000000 -append \\\n")
        if STrun:
            cont_file.write("-s ST_{}.tpr \\\n".format(filename))
            cont_file.write("-cpi ST_{}.cpt \\\n".format(filename))
            cont_file.write("-deffnm ST_{}\n".format(filename))
        elif MDrun:
            cont_file.write("-s md_{}.tpr \\\n".format(filename))
            cont_file.write("-cpi md_{}.cpt \\\n".format(filename))
            cont_file.write("-deffnm md_{}\n".format(filename))
        else:
            pass
        cont_file.close()
        if MDrun:
            repeat_file = open((filename+'MDTitanrepeat.sh'), 'w')
            repeat_file.write('export item=`qsub {}MDTitanpbs.pbs`\n'.format(filename))
            repeat_file.write('for i in {0..15}\n')
            repeat_file.write('do\n')
            repeat_file.write("     item=$(qsub -W depend=afterany:$item {}MDTitancont.pbs) \n".format(filename))
            repeat_file.write('done\n')
        elif STrun:
            repeat_file = open((filename+'STTitanrepeat.sh'), 'w')
            repeat_file.write('export item=`qsub {}STTitanpbs.pbs`\n'.format(filename))
            repeat_file.write('for i in {0..15}\n')
            repeat_file.write('do\n')
            repeat_file.write("     item=$(qsub -W depend=afterany:$item {}STTitancont.pbs) \n".format(filename))
            repeat_file.write('done\n')


    def write_Accre_script(self, STrun = False, MDrun = False):
        filename = self._filename
        if STrun:
            init_file = open((filename + 'STAccresbatch.sbatch'),'w')
        elif MDrun:
            init_file = open((filename + 'MDAccresbatch.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        init_file.write('#!/bin/bash \n')
        init_file.write('#SBATCH --nodes=1\n')
        init_file.write('#SBATCH --account=mccabe_gpu\n')
        init_file.write('#SBATCH --ntasks-per-node=12\n')
        init_file.write('#SBATCH --partition=maxwell\n')
        init_file.write('#SBATCH --gres=gpu:4\n')
        if STrun:
            init_file.write('#SBATCH --time 60:00:00\n')
        elif MDrun:
            init_file.write('#SBATCH --time 60:00:00\n')
        else:
            pass
        init_file.write('#SBATCH --output=my.stdout\n')
        init_file.write('#SBATCH --mail-type=ALL\n')
        init_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        init_file.write('setpkgs -a gromacs_5.1.2_roce\n')
        init_file.write('cd ~/Trajectories/{}/ \n'.format(filename))
        init_file.write('export OMP_NUM_THREADS=1\n')
        if STrun:
            init_file.write('srun -n 12 gmx_mpi mdrun -ntomp 1 -gpu_id 000111222333 -deffnm ST_{} >& out.log\n'.format(filename))
        elif MDrun:
            init_file.write('srun -n 12 gmx_mpi mdrun -ntomp 1 -gpu_id 000111222333 -deffnm md_{} >& out.log\n'.format(filename))
        else:
            pass
            
        init_file.close()
        
        # Write continue script
        if STrun:
            cont_file = open((filename + 'STAccrecont.sbatch'),'w')
        elif MDrun:
            cont_file = open((filename + 'MDAccrecont.sbatch'),'w')
        else:
            sys.exit("Specify ST or MD")
        cont_file.write('#!/bin/bash \n')
        cont_file.write('#SBATCH --nodes=1\n')
        cont_file.write('#SBATCH --account=mccabe_gpu\n')
        cont_file.write('#SBATCH --ntasks-per-node=12\n')
        cont_file.write('#SBATCH --partition=maxwell\n')
        cont_file.write('#SBATCH --gres=gpu:4\n')
        if STrun:
            cont_file.write('#SBATCH --time 60:00:00\n')
        elif MDrun:
            cont_file.write('#SBATCH --time 60:00:00\n')
        else:
            pass
        cont_file.write('#SBATCH --output=my.stdout\n')
        cont_file.write('#SBATCH --mail-type=ALL\n')
        cont_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
        cont_file.write('setpkgs -a gromacs_5.1.2_roce\n')
        cont_file.write('cd ~/Trajectories/{}/ \n'.format(filename))
        cont_file.write('export OMP_NUM_THREADS=1\n')
        if STrun:
            cont_file.write('srun -n 12 gmx_mpi mdrun -ntomp 1 -gpu_id 000111222333 -append -cpi ST_{}.cpt \\\n'.format(filename))
            cont_file.write('-s ST_{}.tpr \\\n'.format(filename))
            cont_file.write('-deffnm ST_{} >& out.log\n'.format(filename))
        elif MDrun:
            cont_file.write('srun -n 12 gmx_mpi mdrun -ntomp 1 -gpu_id 000111222333 -append -cpi md_{}.cpt \\\n'.format(filename))
            cont_file.write('-s md_{}.tpr \\\n'.format(filename))
            cont_file.write('-deffnm md_{} >& out.log\n'.format(filename))
        else:
            pass
        
        cont_file.close()
        
        
        # Write repeat script
        if STrun:
            repeat_file = open((filename + 'STrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}STAccresbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}STAccrecont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}STAccrecont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        elif MDrun:
            repeat_file = open((filename + 'MDrepeat.sh'),'w')
            repeat_file.write('export item=`sbatch {}MDAccresbatch.sbatch` \n'.format(filename))
            repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}MDAccrecont.sbatch`\n'.format(filename))
            repeat_file.write('for i in {0..0}\n')
            repeat_file.write('do\n')
            repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}MDAccrecont.sbatch)\n'.format(filename))
            repeat_file.write('done')
        else:
            pass
        repeat_file.close()

