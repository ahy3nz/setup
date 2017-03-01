import mbuild as mb
import warnings
import pdb
import numpy as np
import mdtraj as mdtraj
import sys
from optparse import OptionParser
from Prototypes_CG import *

def make_layer(n_x = 8, n_y = 8, lipids = None, tilt_angle = 0, spacing = 0, layer_shift = 0, res_index = 0, table_of_contents = None,
        random_z_displacement = 0):
    """ Generate a bilayer leaflet by randomly laying down leaflets in a 2D grid

    Parameters
    ---------
    n_x : int
        2D grid dimension
    n_y : int
        2D grid dimension
    tilt_angle : float
        tilt angle (spun around y-axis)
    spacing : float
        spacing between leaflets, based on area per lipid
    layer_shift : float
        Leaflet weparation from z-axis (used to separate bilayer leaflets)
    res_index : int
        Starting residue index for leaflet construction and residue counting
    table_of_contents : file
        Output file listing residue index, residue name, n_particles for tha residue
    random_z_displacement : float
        Randomly offset molecules by a small amount
    
    Returns
    -------
    layer : mb.Compound()
        Leaflet of molecules
    resindex : int
        Running count of molecules placed (excluding waters)
  
    """
    layer = mb.Compound()
    for i in range(n_x):
        for j in range(n_y):
            # Randomly select a lipid that has not yet been selected
            random_lipid = np.random.randint(0, len(lipids))

            # Create the mbuild Molecule for this lipid
            molecule_i = lipids.pop(random_lipid)
            molecule_to_add = mb.clone(molecule_i[0])

            # Apply tilt angle
            mb.spin_y(molecule_to_add, tilt_angle)

            # Apply z_offset
            z_offset = molecule_i[1]

            # Apply APL and z_offset to identify the position for the molecule in the grid
            position = [i * spacing, j * spacing, z_offset + 
                    (-1 * np.random.random() * random_z_displacement)]
            mb.translate(molecule_to_add, position)

            # Add the new molecule to the layer
            layer.add(molecule_to_add)
            
            # Add to table of contents
            table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, molecule_i[0].name, molecule_i[0].n_particles))
            res_index += 1
    # Apply final layer shifts
    mb.translate(layer, [0, 0, layer_shift])
    #mb.translate(water_layer, [0, 0, layer_shift])
    #return layer, res_index, water_layer
    return layer, res_index


"""
    # Construction of water layer goes underneath the current layer
    water_layer = mb.Compound()
    if (n_solvent_per_lipid):
        for i in range(n_x):
          for j in range(n_y):
              for k in range(n_solvent_per_lipid):
                 water_to_add = water() 
                 # Add waters close to the layer, then build downwards under the layer
                 z_offset = 2
                 water_spacing = 1.5 
                 position = [i * spacing, j * spacing, z_offset - (k * water_spacing) ]
                 mb.translate(water_to_add, position)
                 water_layer.add(water_to_add)
    # Shift all coordinates such that everything lies in positive Z
    min_z_shift = min( [ min(layer.xyz[:,2]), min(water_layer.xyz[:,2])])
    print(min_z_shift)
    mb.translate(layer, [0, 0, -1 * min_z_shift])
    mb.translate(water_layer, [0, 0, -1 * min_z_shift])
"""
    

def enumerate_system_components(lipid_system_info = None):
    """ Generate lipid leaflet definitions
    
    Parameters
    ----------
    lipid_system_info : list
        List whose entries are lists of (mb.Molecule, # of that molecule, molecule z offset)

    Returns
    ------
    top_lipids : list
        list whose entries correspond to each lipid (mb.Molecule(), z_offset)
    bot_lipids : list
        list whose entries correspond to each lipid (mb.Molecule(), z_offset)

    """
    top_lipids = []
    bot_lipids = []

    # Loop through each type of molecule (DSPC, DPPC, etc.)
    for i, lipid_type in enumerate(lipid_system_info):
        n_molecule_per_leaflet = int(lipid_type[1]/2)
        # Loop through the system's quantity of that particular molecule
        #for n in range(int(lipid_type[1])):
        for n in range(n_molecule_per_leaflet):
            lipid_molecule = lipid_type[0]
            lipid_offset = lipid_type[2]
            top_lipids.append([lipid_molecule, lipid_offset])
            bot_lipids.append([lipid_molecule, lipid_offset])
    return top_lipids, bot_lipids


def write_top_file(filename = 'default', lipid_system_info = None):
    """ Generate topology file

    Parameters
    ----------
    filename : Str
        Filename for topology file
    lipid_system_info : list
        List whose entries are lists of (mb.Molecule, # of that molecule, molecule z offset)

    Notes:
    ------
    Assumes a particular path for itp files

    """

    top_file = open(filename + '.top', 'w')

    # Include statment, edit for path to maritini FF
    top_file.write("#include /raid6/homes/ahy3nz/Programs/setup/FF/CG/martini_ff.itp \n")
    top_file.write("\n[ system ]\n")
    top_file.write("Coarse-grained bilayer system\n")
    top_file.write("\n[ molecules ] \n") 

    # Using lipid system information, iterate through each molecule, printing out name and number of molecules 
    for i, lipid_type in enumerate(lipid_system_info):
        n_molecule = int(lipid_type[1])
        molecule_name = lipid_type[0].name
        if n_molecule != 0:
            top_file.write("{:<10s}{:<10d}\n".format(molecule_name, n_molecule))

def solvate_bilayer(system = None, n_x = 8, n_y = 8, n_solvent_per_lipid = 5, water_spacing = 0.8, 
        res_index = 0, table_of_contents = None):
    """ Solvate the top and bottom parts of the bilayer

    Parameters
    ---------
    n_x : int
        3D grid dimension
    n_y : int
        3D grid dimension
    n_solvent_per_lipid : int
        Number of water beads per lipid, another grid dimension
    res_index : int
        Starting residue index for leaflet construction and residue counting
    table_of_contents : file
        Output file listing residue index, residue name, n_particles for tha residue


    Returns
    -------
    solvated_system : mb.Compound()
        System with water solvating the outside of the bilayer

    """
    # Construct two 3D grids of water
    # Add to table of contents file for post processing
    # Compute distances to translate such that water is either below or above bilayer
    cube = mb.Grid3DPattern(n_x, n_y, n_solvent_per_lipid)
    cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    bot_water_list = cube.apply(water())
    for i in range(n_x * n_y * n_solvent_per_lipid):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, water().name, water().n_particles))
        res_index += 1 
    bot_water = mb.Compound()
    for compound in bot_water_list:
        bot_water.add(compound)
    highest_botwater = max(bot_water.xyz[:,2])
    lowest_botlipid = min(system.xyz[:,2])
    shift_botwater = abs(highest_botwater - lowest_botlipid) + 1
    mb.translate(bot_water, [0, 0, -1 * shift_botwater])

    cube = mb.Grid3DPattern(n_x, n_y, n_solvent_per_lipid)
    cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    top_water_list = cube.apply(water())
    for i in range(n_x * n_y * n_solvent_per_lipid):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, water().name, water().n_particles))
        res_index += 1
    top_water = mb.Compound()
    for compound in top_water_list:
        top_water.add(compound)
    lowest_topwater = min(top_water.xyz[:,2])
    highest_toplipid = max(system.xyz[:,2])
    shift_topwater = abs(highest_toplipid - lowest_topwater) + 1
    mb.translate(top_water, [0, 0, shift_topwater])
    pdb.set_trace()
    system.add(bot_water)
    system.add(top_water)

    return system

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "CG_bilayer", dest = "filename")
parser.add_option("-a", "--APL", action="store",type="float", default = 0.36, dest = "area_per_lipid")
parser.add_option("-r", "--rot", action="store", type ="float", default = 0.0, dest = "rotation")
parser.add_option("--DSPC", action="store",type="float", default = 0.5, dest = "DSPC_frac")
parser.add_option("--DPPC", action="store",type="float", default = 0.0, dest = "DPPC_frac")
parser.add_option("--C16FFA", action="store",type="float", default = 0.0, dest = "C16FFA_frac")
parser.add_option("--C22FFA", action="store",type="float", default = 0.0, dest = "C22FFA_frac")
parser.add_option("--C12OH", action="store",type="float", default = 0.5,  dest = "C12OH_frac")
parser.add_option("--C14OH", action="store",type="float", default = 0.0, dest = "C14OH_frac")
parser.add_option("--C16OH", action="store",type="float", default = 0.0, dest = "C16OH_frac")
parser.add_option("--C18OH", action="store",type="float", default = 0.0, dest = "C18OH_frac")
parser.add_option("--C20OH", action="store",type="float", default = 0.0, dest = "C20OH_frac")
parser.add_option("--C22OH", action="store",type="float", default = 0.0, dest = "C22OH_frac")
parser.add_option("--C24OH", action="store",type="float", default = 0.0, dest = "C24OH_frac")
parser.add_option("--ISIS", action="store",type="float", default = 0.0, dest = "ISIS_frac")
parser.add_option("--SS", action="store",type="float", default = 0.0, dest = "SS_frac")
parser.add_option("--CHOL", action="store",type="float", default = 0.0, dest = "CHOL_frac")
parser.add_option("--PMEA", action="store",type="float", default = 0.0, dest = "PMEA_frac")
#parser.add_option("--water", action="store",type="float", default = 0.0, dest = "Water_frac")
(options, args) = parser.parse_args()


filename = options.filename
tilt_angle = options.rotation * np.pi/180
area_per_lipid = options.area_per_lipid
spacing = np.sqrt(area_per_lipid)

# Default parameters, less likely to be changed
n_x = 8
n_y = 8
n_lipid = 2 * n_x * n_y
n_solvent_per_lipid = 5 # This is usually 20 waters per molecule, but a water bead is 4 waters
n_solvent = n_lipid * n_solvent_per_lipid
random_z_displacement = 0.3


"""
# Todo: Obtain prototypes for more bilayer molecules
lipid_system_info = [(DSPC(), np.ceil(n_lipid * options.DSPC_frac), 3.2),
                      (DPPC(), np.ceil(n_lipid * options.DPPC_frac), 2.6)
                      (C16FFA(), np.floor(n_lipid * options.C16FFA_frac), 4.4),
                      (C22FFA(), np.floor(n_lipid * options.C22FFA_frac), 3.2),
                      (C12OH(), np.floor(n_lipid * options.C12OH_frac), 2.6),
                      (C14OH(), np.floor(n_lipid * options.C14OH_frac), 2.8),
                      (C16OH(), np.floor(n_lipid * options.C16OH_frac), 2.6),
                      (C18OH(), np.floor(n_lipid * options.C18OH_frac), 2.4),
                      (C20OH(), np.floor(n_lipid * options.C20OH_frac), 3.0),
                      (C22OH(), np.floor(n_lipid * options.C22OH_frac), 3.2),
                      (C24OH(), np.floor(n_lipid * options.C24OH_frac), 3.2),
                      (ISIS(), np.floor(n_lipid * options.ISIS_frac), 3.0),
                      (CHOL(), np.floor(n_lipid * options.CHOL_frac), 4.0)] 
                      """

lipid_system_info = [(DSPC(), np.ceil(n_lipid * options.DSPC_frac), 3.2),
                     (C12OH(), np.floor(n_lipid * options.C12OH_frac), 2.6)]
                     #(water(), np.floor(n_lipid * n_solvent_per_lipid), 6.0)]
                
# Create a list whose entries correspond to a single lipid in the system
top_lipids, bot_lipids = enumerate_system_components(lipid_system_info = lipid_system_info)

if len(bot_lipids)!= (n_lipid)/2:
    sys.exit('Error setting up system components')

# Table of contents to externally store residue indices and names for post setup
table_of_contents = open(filename+'.dat', 'w')
res_index = 1

# Generate bottom layer randomly
bot_layer, res_index = make_layer(n_x = 8, n_y = 8, lipids = bot_lipids, tilt_angle = tilt_angle, spacing = spacing, layer_shift = 0,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement)

# Generate the top layer randomly
top_layer, res_index = make_layer(n_x = 8, n_y = 8, lipids = top_lipids, tilt_angle = tilt_angle, spacing = spacing, layer_shift = 3.2,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement)
       
# Rotate bottom layer to form bilayer
mb.spin_y(bot_layer, theta=np.pi)

# Create system class that includes top and bottom layers
system = mb.Compound()
system.add(bot_layer)
system.add(top_layer)

# Solvate system
system = solvate_bilayer(system = system, n_x = n_x, n_y = n_y, n_solvent_per_lipid = n_solvent_per_lipid, 
        res_index = res_index, table_of_contents = table_of_contents)


# Shift everything to positive z
min_z_shift = min(system.xyz[:,2])
mb.translate(system, [0, 0, -1 * min_z_shift])

# Write gro file suppressing mbuild warnings
print("Writing <{0}.gro> ...".format(filename))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    system.save(filename + '.gro',overwrite=True)
table_of_contents.close()

# Write topology file
print("Writing <{0}.top> ...".format(filename))
write_top_file(filename = filename, lipid_system_info = lipid_system_info)

