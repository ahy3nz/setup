import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")

(options, args) = parser.parse_args()
filename = options.filename

# Write initial job script
init_file = open((filename + 'sbatch.sbatch'),'w')
init_file.write('#!/bin/bash -l\n')
init_file.write('#SBATCH -p regular\n')
init_file.write('#SBATCH -N 8\n')
init_file.write('#SBATCH -t 10:00:00\n')
init_file.write('#SBATCH -J {}\n'.format(filename))
init_file.write('#SBATCH -o {}\n'.format(filename))
init_file.write('#SBATCH -C haswell\n')
init_file.write('#SBATCH --mail-type=ALL\n')
init_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
init_file.write('module load gromacs/5.1.2\n')
init_file.write('srun -n 128 -c 4 mdrun_mpi_sp -deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))

# Write continue script
cont_file = open((filename + 'cont.sbatch'),'w')
cont_file.write('#!/bin/bash -l\n')
cont_file.write('#SBATCH -p regular\n')
cont_file.write('#SBATCH -N 8\n')
cont_file.write('#SBATCH -t 10:00:00\n')
cont_file.write('#SBATCH -J {}\n'.format(filename))
cont_file.write('#SBATCH -o {}\n'.format(filename))
cont_file.write('#SBATCH --mail-type=ALL\n')
cont_file.write('#SBATCH -C haswell\n')
cont_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
cont_file.write('module load gromacs/5.1.2\n')
cont_file.write('srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/{}/md_{}.cpt \\\n'.format(filename,filename))
cont_file.write('-s $SCRATCH/Trajectories/{}/md_{}.tpr \\\n'.format(filename,filename))
cont_file.write('-deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
cont_file.close()




# Write repeat script
repeat_file = open((filename + 'repeat.sh'),'w')
repeat_file.write('export item=`sbatch {}sbatch.sbatch` \n'.format(filename))
repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}cont.sbatch`\n'.format(filename))
repeat_file.write('for i in {0..0}\n')
repeat_file.write('do\n')
repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}cont.sbatch)\n'.format(filename))
repeat_file.write('done')
repeat_file.close()
