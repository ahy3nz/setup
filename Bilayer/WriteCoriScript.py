import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")
parser.add_option("--ST", action="store_true", dest = "STrun")
parser.add_option("--MD", action="store_true", dest = "MDrun")

(options, args) = parser.parse_args()
filename = options.filename

# Write initial job script
if options.STrun:
    init_file = open((filename + 'STsbatch.sbatch'),'w')
elif options.MDrun:
    init_file = open((filename + 'MDsbatch.sbatch'),'w')
else:
    sys.exit("Specify ST or MD")
init_file.write('#!/bin/bash -l\n')
init_file.write('#SBATCH -p regular\n')
init_file.write('#SBATCH -N 8\n')
init_file.write('#SBATCH -t 09:00:00\n')
init_file.write('#SBATCH -J {}\n'.format(filename))
init_file.write('#SBATCH -o {}\n'.format(filename))
init_file.write('#SBATCH -C haswell\n')
init_file.write('#SBATCH --mail-type=ALL\n')
init_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
init_file.write('module load gromacs/5.1.2\n')
if options.STrun:
    init_file.write('srun -n 128 -c 4 mdrun_mpi_sp -deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
elif options.MDrun:
    init_file.write('srun -n 128 -c 4 mdrun_mpi_sp -deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
else:
    pass
    
init_file.close()

# Write continue script
if options.STrun:
    cont_file = open((filename + 'STcont.sbatch'),'w')
elif options.MDrun:
    cont_file = open((filename + 'MDcont.sbatch'),'w')
else:
    sys.exit("Specify ST or MD")
cont_file.write('#!/bin/bash -l\n')
cont_file.write('#SBATCH -p regular\n')
cont_file.write('#SBATCH -N 8\n')
cont_file.write('#SBATCH -t 09:00:00\n')
cont_file.write('#SBATCH -J {}\n'.format(filename))
cont_file.write('#SBATCH -o {}\n'.format(filename))
cont_file.write('#SBATCH --mail-type=ALL\n')
cont_file.write('#SBATCH -C haswell\n')
cont_file.write('#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu\n')
cont_file.write('module load gromacs/5.1.2\n')
if options.STrun:
    cont_file.write('srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/{}/ST_{}.cpt \\\n'.format(filename,filename))
    cont_file.write('-s $SCRATCH/Trajectories/{}/ST_{}.tpr \\\n'.format(filename,filename))
    cont_file.write('-deffnm $SCRATCH/Trajectories/{}/ST_{} >& out.log\n'.format(filename,filename))
elif options.MDrun:
    cont_file.write('srun -n 128 -c 4 mdrun_mpi_sp -append -cpi $SCRATCH/Trajectories/{}/md_{}.cpt \\\n'.format(filename,filename))
    cont_file.write('-s $SCRATCH/Trajectories/{}/md_{}.tpr \\\n'.format(filename,filename))
    cont_file.write('-deffnm $SCRATCH/Trajectories/{}/md_{} >& out.log\n'.format(filename,filename))
else:
    pass

cont_file.close()


# Write repeat script
if options.STrun:
    repeat_file = open((filename + 'STrepeat.sh'),'w')
    repeat_file.write('export item=`sbatch {}STsbatch.sbatch` \n'.format(filename))
    repeat_file.write('#export item=`sbatch --dependency=afterany:2418639 {}STcont.sbatch`\n'.format(filename))
    repeat_file.write('for i in {0..0}\n')
    repeat_file.write('do\n')
    repeat_file.write('	item=$(sbatch --dependency=afterany:${{item:20:7}} {}STcont.sbatch)\n'.format(filename))
    repeat_file.write('done')
elif options.MDrun:
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
