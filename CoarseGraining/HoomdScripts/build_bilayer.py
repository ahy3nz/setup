import numpy as np

import mbuild as mb
from mbuild.formats.hoomdxml import write_hoomdxml
from foyer import Forcefield

import scripts.bilayer as bilayer

# Import statements for molecule prototypes
#import atomistic.ecer2_hairpin.ecer2 as ecer2
#import atomistic.c24ffa.ffa24 as ffa24
#import atomistic.tip3p.SOL as SOL
from cg.prototypes.DSPC import DSPC

###################
## Sample script to construct and save an mBuild Bilayer for Gromacs
####################

# Bilayer specifications
filename = 'compound'
HOOMD_FF = '/raid6/homes/ahy3nz/Programs/setup/FF/CG/msibi_ff.xml'
apl = 0.54
tilt_angle = np.deg2rad(5)
solvent_density = 900
random_spin=np.deg2rad(20)
n_x = 8
n_y = 8
n_solvent_per_lipid = 10
leaflet_info = [ (DSPC(), 64, 0)
                 ]

system = bilayer.Bilayer(leaflet_info=leaflet_info, n_x=n_x, n_y=n_y, apl=apl, 
        tilt_angle=tilt_angle, random_spin=random_spin,
        solvent=mb.Particle(name="_W"), solvent_density=solvent_density,
        n_solvent_per_lipid=n_solvent_per_lipid, solvent_mass=72)

#system = bilayer.translate_to_positive_octant(system)
system = bilayer.center_on_origin(system)
box = mb.Box(lengths=system.boundingbox.lengths)
structure = system.to_parmed(box=box)

ff = Forcefield(forcefield_files=HOOMD_FF)
kwargs  = {}
kwargs['rigid_bodies'] = [p.rigid_id for p in system.particles()]
structure = ff.apply(structure, assert_dihedral_params=False)
write_hoomdxml(structure, '{}.hoomdxml'.format(filename), 
        ref_energy = 0.239, ref_distance = 1, **kwargs)

