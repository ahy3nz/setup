import hoomd
import hoomd.md
import numpy as np
import pdb
import itertools
from optparse import OptionParser
import martini_hoomd

import hoomd.deprecated
hoomd.context.initialize("");

# load a structure, generated from mbuild
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "CG_bilayer", dest = "filename")
(options, args) = parser.parse_args()
filename = options.filename
# load a structure, generated from mbuild
#hoomd.init.read_gsd(filename, restart= 'trajectory.gsd')
hoomd.init.read_gsd(filename, restart = 'restart.gsd')


atomtypes = martini_hoomd.get_atom_types()
hoomd.util.quiet_status()
martini_hoomd.set_bonds()
#martini_hoomd.set_harmonic_angles(timtails = True) # For using harmonic-fit potentials
martini_hoomd.set_cosine_angles(timtails = True)
nl = hoomd.md.nlist.cell()

# Set exclusions only if using martinitim
martini_hoomd.set_exclusions(nl=nl)
lamb0 = martini_hoomd.set_pairs(table_dir = "lambda0.9-martinitim",nl=nl)
#lamb0 = martini_hoomd.set_pairs(table_dir = "lambda0.9-martinifull",nl=nl)

# Define some groups
all = hoomd.group.all()
hoomd.util.unquiet_status()


# Output parameters
hoomd.analyze.log(filename = "em1.log", 
        quantities = ['potential_energy', 'temperature', 'lx', 'ly', 'lz'],
    period=100, overwrite=True);
# If this is a continuation script, just append to the trajectory
hoomd.dump.gsd("trajectory.gsd", period=500, group=all, phase=0);
# Restart files should overwrite and replace frame 0
hoomd.dump.gsd("restart.gsd", period=500, group=all, overwrite=True, truncate=True,phase=0);



# FIRE algorithm for EM
hoomd.md.integrate.mode_standard(dt=0.005);
EM = hoomd.md.integrate.mode_minimize_fire(group = all, dt = 0.005)
hoomd.run(10000); # number of steps in time
