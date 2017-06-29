import hoomd
import hoomd.md
import martini_hoomd
from optparse import OptionParser

import hoomd.deprecated
hoomd.context.initialize("");

# load a structure, generated from mbuild
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "CG_bilayer", dest = "filename")
(options, args) = parser.parse_args()
filename = options.filename
# load a structure, generated from mbuild
hoomd.init.read_gsd(filename, restart= 'restart.gsd')
#hoomd.init.read_gsd('trajectory.gsd')

atomtypes = martini_hoomd.get_atom_types()
ENERGY_TABLE, DISTANCE_TABLE, FORCE_TABLE = martini_hoomd.get_table_units()
hoomd.util.quiet_status()
martini_hoomd.set_bonds()
#martini_hoomd.set_harmonic_angles(timtails = True) # For using harmonic-fit potentials
martini_hoomd.set_cosine_angles(timtails = True)
nl = hoomd.md.nlist.cell()

# Set exclusions only if using martinitim
martini_hoomd.set_exclusions(nl=nl)
lamb0 = martini_hoomd.set_pairs(table_dir = "lambda0-martinitim",nl=nl)
#lamb0 = martini_hoomd.set_pairs(table_dir = "lambda0-martinifull",nl=nl)

# Set reaction field parameters
#martini_hoomd.set_reactionfield_electrostatics(nl=nl)

# Define some groups
all = hoomd.group.all()
hoomd.util.unquiet_status()


# Output parameters
hoomd.analyze.log(filename = "nvtnpt.log", 
        quantities = ['potential_energy', 'temperature', 'lx', 'ly', 'lz'],
    period=100, overwrite=True);

# If this is a continuation script, just append to the trajectory
hoomd.dump.gsd("trajectory.gsd", period=500, group=all, phase=0);
# Restart files should overwrite and replace frame 0
hoomd.dump.gsd("restart.gsd", period=500, group=all, overwrite=True, truncate=True,phase=0);
hoomd.deprecated.dump.xml(all, filename="restart.hoomdxml", period=500, restart=True,phase=0);




# NVT
#hoomd.md.integrate.mode_standard(dt=0.005);
#nvtintegrator = hoomd.md.integrate.nvt(group=all, kT=2.5, tau = 0.2);
#hoomd.run(10000); # number of steps in time

# NPT cooling 
total_timesteps = 1e8
cooling_variant = hoomd.variant.linear_interp(points = [(0,3), (1e8, 1)])
hoomd.md.integrate.mode_standard(dt=0.01)
cooling_integrator = hoomd.md.integrate.npt(group = all, couple="xy", kT=cooling_variant, tau = 0.2, P = 0.060221, tauP = 0.5)
nptsnapshot = hoomd.dump.gsd("npt.gsd", period=500, group=all, overwrite=True, truncate=True, phase=0);
hoomd.run(total_timesteps)
cooling_integrator.disable()

# NPT heating
total_timesteps = 1e8
heating_variant = hoomd.variant.linear_interp(points = [(0,1), (1e8, 3)])
hoomd.md.integrate.mode_standard(dt=0.01)
heating_integrator = hoomd.md.integrate.npt(group = all, couple="xy", kT=heating_variant, tau = 0.2, P = 0.060221, tauP = 0.5)
hoomd.run(total_timesteps)

