import hoomd
import hoomd.md
import martini_hoomd
from optparse import OptionParser

import hoomd.deprecated
hoomd.context.initialize("");
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default="CG_bilayer", dest="filename")
parser.add_option("-t", action="store_true", dest="timtails", default=False)
(options, args) = parser.parse_args()
filename = options.filename
# load a structure, generated from mbuild
hoomd.init.read_gsd(filename, restart= 'restart.gsd')
#hoomd.init.read_gsd('trajectory.gsd')


atomtypes = martini_hoomd.get_atom_types()
hoomd.util.quiet_status()
martini_hoomd.set_bonds()
#martini_hoomd.set_harmonic_angles(timtails = True) # For using harmonic-fit potentials
martini_hoomd.set_cosine_angles(timtails=options.timtails)
nl = hoomd.md.nlist.cell()

# Set exclusions only if using martinitim
martini_hoomd.set_exclusions(nl=nl)
if options.timtails:
    table_dir = "lambda0.25-martinitim"
else:
    table_dir = "lambda0.25-martinifull"
lamb0 = martini_hoomd.set_pairs(table_dir=table_dir, nl=nl)


# Set reaction field parameters
#martini_hoomd.set_reactionfield_electrostatics(nl=nl)

# Define some groups
all = hoomd.group.all()
hoomd.util.unquiet_status()

# Output parameters
hoomd.analyze.log(filename = "em3.log", 
        quantities = ['potential_energy', 'temperature', 'lx', 'ly', 'lz'],
    period=500, overwrite=True);

# If this is a continuation script, just append to the trajectory
hoomd.dump.gsd("trajectory.gsd", period=500, group=all, phase=0);
# Restart files should overwrite and replace frame 0
hoomd.dump.gsd("restart.gsd", period=500, group=all, overwrite=True, truncate=True,phase=0);


# FIRE algorithm for EM, now including charges
hoomd.md.integrate.mode_standard(dt=0.005);
EM = hoomd.md.integrate.mode_minimize_fire(group = all, dt = 0.005)
hoomd.run(10000); # number of steps in time
