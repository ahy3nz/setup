import os
import numpy as np
import hoomd
import hoomd.md
from hoomd.deprecated.init import read_xml
import hoomd_initialization
import itertools

hoomd.context.initialize("")
system = read_xml(filename="compound.hoomdxml", wrap_coordinates=True)
kb = 8.314e-3
atm = 0.060221

pot_width = 121
nl = hoomd.md.nlist.cell()
table = hoomd.md.pair.table(width=pot_width, nlist=nl)
potential_dir = '/scratch/yangah/MSIBI_FF/FF1.17.28/potentials/'

beadtypes = ['W', 'E1','C2', 'C3', 'PCP', 'PCN']
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
hoomd.util.unquiet_status()
zero_momentum = hoomd.md.update.zero_momentum(period=1e4)
hoomd.analyze.log('thermo.log', ['lx', 'ly', 'lz', 'temperature', 'pressure', 
    'potential_energy'], overwrite=True, period=1000)
dcd_traj = hoomd.dump.dcd('query.dcd', period=10000, group=all, 
        phase=0, overwrite=True, 
        unwrap_full=False)
gsd_restart = hoomd.dump.gsd(filename='restart.gsd', group=all, 
        truncate=True, period=1e6)

# NVT with changing box
box_dims = [system.box.get_lattice_vector(0)[0],
        system.box.get_lattice_vector(1)[1],
        system.box.get_lattice_vector(2)[2]]
lx_variant = hoomd.variant.linear_interp([
    (0, box_dims[0]), (20e6, box_dims[0] * np.sqrt(2)), 
    (40e6, box_dims[0])
    ])
ly_variant = hoomd.variant.linear_interp([
    (0, box_dims[1]), (20e6, box_dims[1] * np.sqrt(2)), 
    (40e6, box_dims[1])
    ])
lz_variant = hoomd.variant.linear_interp([
    (0, box_dims[2]), (20e6, box_dims[2] / 2), 
    (40e6, box_dims[2])
    ])

temp_variant = hoomd.variant.linear_interp([
    (0, 500*kb), (20e6, 500*kb),
    (40e6, 405*kb)
    ])
hoomd.update.box_resize(Lx=lx_variant, Ly=ly_variant, Lz=lz_variant)

hoomd.md.integrate.mode_standard(dt=0.02)
nvt_integrator = hoomd.md.integrate.nvt(group=all, kT=temp_variant, tau=1)
hoomd.run_upto(40e6)
nvt_integrator.disable()

# NPT with annealing
temp_variant = hoomd.variant.linear_interp([
    (40e6, 450*kb),
    (60e6, 305*kb)
    ])
hoomd.md.integrate.mode_standard(dt=0.02)
npt_integrator = hoomd.md.integrate.npt(group=all, couple='xy', 
        kT=temp_variant, tau=1,
        P=atm, tauP=10)
hoomd.run_upto(60e6)
gsd_restart.write_restart()


