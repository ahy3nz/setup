import mbuild as mb
import pdb
import numpy as np
import mdtraj as mdtraj
import sys
from DSPC_CG import *
from C12OH_CG import *


n_x = 8
n_y = 8
n_lipid = n_x * n_y
n_solvent_per_layer = 1280
DSPC_frac = 1.0
C12OH_frac = 0.0
DPPC_frac, C14OH_frac, C16OH_frac, C18OH_frac, C20OH_frac, C22OH_frac, C24OH_frac = 0,0,0,0,0,0,0
C22FFA_frac, C16FFA_frac, ISIS_frac, SS_frac, CHOL_frac, PMEA_frac = 0,0,0,0,0,0
tilt_angle = 0*np.pi/180
area_per_lipid = 0.3
spacing=np.sqrt(area_per_lipid)


# I could use mBuild to construct the system with their geometric transformations, but then i'd need to atomtype it myself?
# Could I just create a template of the entire bilayer, with correct atomtypes, and then fill in the coordinates from mBuild?
"""
# Or would it be easier to use mBuild to build bilayer system, and update/replace that with my own atomtypes ( probably this)

lipid_system_info = [(DSPC_CG(), np.ceil(n_lipid * DSPC_frac), 3.2),
                      (DPPC_CG(), np.ceil(n_lipid * DPPC_frac), 2.6)
                      (C16FFA(), np.floor(n_lipid * C16FFA_frac), 4.4),
                      (C22FFA(), np.floor(n_lipid * C22FFA_frac), 3.2),
                      (C12OH(), np.floor(n_lipid * C12OH_frac), 2.6),
                      (C14OH(), np.floor(n_lipid * C14OH_frac), 2.8),
                      (C16OH(), np.floor(n_lipid * C16OH_frac), 2.6),
                      (C18OH(), np.floor(n_lipid * C18OH_frac), 2.4),
                      (C20OH(), np.floor(n_lipid * C20OH_frac), 3.0),
                      (C22OH(), np.floor(n_lipid * C22OH_frac), 3.2),
                      (C24OH(), np.floor(n_lipid * C24OH_frac), 3.2),
                      (ISIS(), np.floor(n_lipid * ISIS_frac), 3.0),
                      (CHOL(), np.floor(n_lipid * CHOL_frac), 4.0)] 
                      """

lipid_system_info = [(DSPC_CG(), np.ceil(n_lipid * DSPC_frac), 3.2),
                     (C12OH_CG(), np.floor(n_lipid * C12OH_frac), 2.6)]
                
# Create a list whose entries correspond to a single lipid in the system
top_lipids = []
bot_lipids = []
index = 0
# Loop through each type of molecule (DSPC, DPPC, etc.)
for i, lipid_type in enumerate(lipid_system_info):
    # Loop through the quantity of that particular molecule
    for n in range(int(lipid_type[1])):
        lipid_molecule = lipid_type[0]
        lipid_offset = lipid_type[2]
        top_lipids.append([lipid_molecule, lipid_offset])
        bot_lipids.append([lipid_molecule, lipid_offset])
        index += 1

if len(bot_lipids)!=n_lipid:
    sys.exit('Error setting up system components')

# Generate bottom layer randomly
bot_layer = mb.Compound()
for i in range(n_x):
    for j in range(n_y):
        # Randomly select a lipid that has not yet been selected
        random_lipid = np.random.randint(0, len(bot_lipids))
        # Create the mbuild Molecule for this lipid
        molecule_i = bot_lipids.pop(random_lipid)
        molecule_to_add = mb.clone(molecule_i[0])
        # Apply tilt angle
        mb.spin_y(molecule_to_add, tilt_angle)
        # Apply z_offset
        z_offset = molecule_i[1]
        #z_offset = 0
        # Apply APL and z_offset to identify the position for the molecule in the grid
        position = [i * spacing, j * spacing, z_offset]
        mb.translate(molecule_to_add, position)
        # Add the new molecule to the layer
        bot_layer.add(molecule_to_add)
       
        
# Generate the top layer of lipids randomly
top_layer = mb.Compound()
for i in range(n_x):
    for j in range(n_y):
        # Randomly select a lipid that has not yet been selected
        random_lipid = np.random.randint(0, len(top_lipids))
        # Create the mbuild Molecule for this lipid
        molecule_i = top_lipids.pop(random_lipid)
        molecule_to_add = mb.clone(molecule_i[0])
        # Apply tilt angle
        mb.spin_y(molecule_to_add, tilt_angle)
        # Apply z_offset
        z_offset = molecule_i[1]
        # Apply APL and z_offset to identify the position for the molecule in the grid
        position = [i * spacing, j * spacing, z_offset+3.2]
        mb.translate(molecule_to_add, position)
        # Add the new molecule to the layer
        top_layer.add(molecule_to_add)
        
# Rotate bottom layer to form bilayer
mb.spin_y(bot_layer, theta=np.pi)
# Create system class that includes top and bottom layers
system = mb.Compound()
system.add(top_layer)
system.add(bot_layer)
system.save('CG_bilayer.gro',overwrite=True)



