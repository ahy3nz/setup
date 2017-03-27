import mbuild as mb
import numpy as np
import pdb
from Prototypes import *


system = ISIS()
mb.translate(system, -system[0].pos)

# Make sure the tilt angle is zero
mb.z_axis_transform(system, point_on_z_axis=system.children[1], point_on_zx_plane=system.children[9])

# If this thing is pointing in the negative z-direction, invert the z-coordiantes so it points positive z
for i in range(len(system.children)):
    system.children[i].pos[2] = np.abs(system.children[i].pos[2])
    system.children[i].pos[1] = np.abs(system.children[i].pos[1])
    system.children[i].pos[0] = np.abs(system.children[i].pos[0])

mb.z_axis_transform(system, point_on_z_axis=system.children[1], point_on_zx_plane=system.children[3])
#mb.spin_x(system, theta=(3.14)*(3/8)) #guessin between 1/2 and 1/4
mb.spin_x(system, theta=(3/8)*3.14)

system.save('Prototypes/new_ISIS.gro',overwrite=True)
