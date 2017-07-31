import os
import subprocess
import pdb
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from optparse import OptionParser
"""
Convert Martini (or any gromacs) LJ parameters 
to tabulated form for HOOMD simulation.
table_energy [=] 0.1 kcal/mol
table_dist [=] 6 angstrom
table_force [=] table_energy/table_dist
specify f, c6, c12
Martini cutoff at 1.1nm
"""
# These are conversion factors to convert from gromacs 
# units (kJ mol-1, nm, kj mol-1  nm-1) to table units
ENERGY_GMX_TABLE = 0.4184
DISTANCE_GMX_TABLE = 0.6
FORCE_GMX_TABLE = ENERGY_GMX_TABLE/DISTANCE_GMX_TABLE

# I just need to multiply the respective entires by the factor

# Specify an input file
parser = OptionParser()
parser.add_option("-f", action = "store", type = "string", dest = "inp")
(options, args) = parser.parse_args()
p = subprocess.Popen('cp {}.txt {}_original.txt'.format(options.inp[ : -4], options.inp[ : -4]), shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
p.wait()
print(p.stdout.read())
tabulated_potential = np.loadtxt(options.inp)
# Scale energy and force terms
tabulated_potential[:,0] *= DISTANCE_GMX_TABLE
tabulated_potential[:,1] *= ENERGY_GMX_TABLE
tabulated_potential[:,2] *= FORCE_GMX_TABLE
# Print to output directory with same filename
#output_file = os.path.join(output_dir, pair_potential)
np.savetxt(options.inp, tabulated_potential)

