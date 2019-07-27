import os
import numpy as np
import hoomd
import hoomd.md
from hoomd.deprecated.init import read_xml
import hoomd_initialization
import itertools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", type=str)
args = parser.parse_args()
config = args.c

hoomd.context.initialize("")
path_to_configs = '/scratch/yangah/MSIBI_FF/FF5.1.15/configs/'
system = read_xml(filename=path_to_configs + config + '.hoomdxml',
        wrap_coordinates=True)
kb = 8.314e-3
atm = 0.060221

pot_width = 121
nl = hoomd.md.nlist.cell()
table = hoomd.md.pair.table(width=pot_width, nlist=nl)
potential_dir = '../potentials/'

beadtypes = ['W', 'HS', 'E1','C2', 'C3', 'PCP', 'PCN']
for i,j in itertools.combinations_with_replacement(beadtypes,2):
    if os.path.isfile('{}/{}-{}.txt'.format(potential_dir,i,j)):
        table.set_from_file(i,j, filename='{}/{}-{}.txt'.format(potential_dir, i,j))
    else:
        table.set_from_file(i,j, filename='{}/{}-{}.txt'.format(potential_dir, j,i))

hoomd.util.quiet_status()
hoomd_initialization.set_bonds(system)
hoomd_initialization.set_harmonic_angles()
hoomd_initialization.set_exclusions(nl=nl)

all = hoomd.group.all()
waters = hoomd.group.type("W", name="W")
nonwaters = hoomd.group.difference("nonW", a=all, b=waters)

hoomd.util.unquiet_status()
zero_momentum = hoomd.md.update.zero_momentum(period=1e4)
dcd_traj = hoomd.dump.dcd('{}.dcd'.format(config), 
        period=10000, group=all, 
        phase=0, overwrite=True, 
        unwrap_full=False)
hoomd.analyze.log('{}.log'.format(config), 
        ['lx', 'ly', 'lz', 
        'potential_energy', 
        'temperature', 'pressure'],
        overwrite=True, period=10000)


# Quick energy minimization
int_mode = hoomd.md.integrate.mode_standard(dt=0.01)
nve_int = hoomd.md.integrate.nve(all, limit=0.5)
print("Running NVE ...")
hoomd.run(4e3)
nve_int.disable()

# NVT thermo annealing, 10fs timestep, 20ns total
t_variant = hoomd.variant.linear_interp([
    (0, 100*kb),
    (5e5, 500*kb),
    (10e5, 500*kb),
    (15e5, 305*kb)
    ], zero='now')
w_nvt_int = hoomd.md.integrate.nvt(group=all, kT=t_variant, tau=1)


print("Running NVT thermo annealing ...")
hoomd.run(20e5)

# NVT box annealing, 10fs timestep, 20ns total
box_dims = [system.box.get_lattice_vector(0)[0],
        system.box.get_lattice_vector(1)[1],
        system.box.get_lattice_vector(2)[2]]
lx_variant = hoomd.variant.linear_interp([
    (0, box_dims[0]), (10e5, box_dims[0] * np.sqrt(2)), 
    (20e5, box_dims[0])
    ], zero='now')
ly_variant = hoomd.variant.linear_interp([
    (0, box_dims[1]), (10e5, box_dims[1] * np.sqrt(2)), 
    (20e5, box_dims[1])
    ], zero='now')
lz_variant = hoomd.variant.linear_interp([
    (0, box_dims[2]), (10e5, box_dims[2] / 2), 
    (20e5, box_dims[2])
    ], zero='now')

box_resizer = hoomd.update.box_resize(Lx=lx_variant, Ly=ly_variant, 
        Lz=lz_variant, period=1e5)
print("Running NVT box annealing ...")
#hoomd.run(20e6)
hoomd.run(5e6)
w_nvt_int.disable()
box_resizer.disable()

# NPT thermo annealing, 20fs timestep
#t_variant = hoomd.variant.linear_interp([
#    (0, 305*kb),
#    (15e5, 305*kb),
#    (30e5, 500*kb),
#    (45e5, 500*kb),
#    (60e5, 305*kb)
#    ], zero='now')
#npt_int = hoomd.md.integrate.npt(group=all, tau=1, kT=t_variant, 
#        P=atm, tauP=10, couple='xy')
#hoomd.run(60e5)
#npt_int.disable()


# NPT equilibation - final
# dt = 0.02
# tau = 1
# tauP = 10
npt_int = hoomd.md.integrate.npt(group=all, tau=1, kT=305*kb, 
        P=atm, tauP=10, couple='xy')
print("Running final NPT ... ")
#hoomd.run(20e6) # 200ns
hoomd.run(5e6) 
