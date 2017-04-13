import mbuild as mb
from collections import OrderedDict
import warnings
import pdb
import numpy as np
import mdtraj as mdtraj
import sys
from itertools import product
from optparse import OptionParser
from Prototypes import *
from scriptWriter import *

GMX_FF_DIR = "/raid6/homes/ahy3nz/Programs/setup/FF/gromos53a6/"

def new_make_layer(n_x = 8, n_y = 8, lipid_system_info = None, tilt_angle = 0, spacing = 0, 
        layer_shift = 0, res_index = 0, table_of_contents = None,
        random_z_displacement = 0, top_file = None, lipid_atom_dict = None, atom_index = 0):
    """ Generate a bilayer leaflet by laying down molecules in a 2D grid at random grid points

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
    lipid_atom_dict : OrderedDict()
        Dictionary whose values are mb.Compounds()s and values are a list of atom indices of that compound
    atom_index : int
        Counter for indexing atoms for lipid_atom_dict
    
    Returns
    -------
    layer : mb.Compound()
        Leaflet of molecules
    resindex : int
        Running count of molecules placed (excluding waters)
    lipid_atom_dict : OrderedDict()
        Dictionary whose values are mb.Compounds()s and values are a list of atom indices of that compound
    atom_index : int
        Counter for indexing atoms for lipid_atom_dict
  
    """


    layer = mb.Compound()

    # Create ordered pairs
    ordered_pairs = []
    for i, j in product(range(n_x), range(n_y)):
        ordered_pairs.append((i,j))


    # Randomly assign ordered pairs to each lipid
    # based on the way lipids is set, all of one molecule is listed first
    # before getting to the next one
    # Loop through each type of molecule (DSPC, DPPC, etc.)
    for i, lipid_type in enumerate(lipid_system_info):
        n_molecule_per_leaflet = int(lipid_type[1]/2)
        #Add to top file
        if(n_molecule_per_leaflet !=0):
            top_file.write("{:<10s}{:<10d}\n".format(lipid_type[0].name, n_molecule_per_leaflet))
        # Loop through the system's quantity of that particular molecule
        for n in range(n_molecule_per_leaflet):
            random_index = np.random.randint(0, len(ordered_pairs))
            (i, j) = ordered_pairs.pop(random_index)

            # Do geometry transformations
            molecule_to_add = mb.clone(lipid_type[0])
            # Apply tilt angle
            mb.spin_y(molecule_to_add, tilt_angle)

            # Apply z_offset
            z_offset = lipid_type[2]

            # Apply APL and z_offset to identify the position for the molecule in the grid
            position = [i * spacing, j * spacing, z_offset + layer_shift +
                        (-1 * np.random.random() * random_z_displacement)]
            mb.translate(molecule_to_add, position)

            # Add the new molecule to the layer
            layer.add(molecule_to_add)
                
            # Add to table of contents
            table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, molecule_to_add.name, molecule_to_add.n_particles))

            # Add to lipid dictionary
            if molecule_to_add.name in lipid_atom_dict:
                lipid_atom_dict[molecule_to_add.name] +=list(range(atom_index, atom_index + molecule_to_add.n_particles, 1))
            else:
                lipid_atom_dict[molecule_to_add.name] = list(range(atom_index, atom_index + molecule_to_add.n_particles,1))

            # Increment counters
            res_index += 1
            atom_index += molecule_to_add.n_particles
        


    return layer, res_index, lipid_atom_dict, atom_index



def write_top_file_header(filename = 'default', lipid_system_info = None):
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

    top_file = open(filename + '2.top', 'w')

    # Include statment, edit for path to maritini FF, make "_b" itp for no charges
    #top_file.write("#include \"{}martini_ff.itp\" \n".format(GMX_FF_DIR))
    #top_file.write(";#include \"{}martini_ff_b.itp\" \n".format(GMX_FF_DIR))
    top_file.write(";#include \"{}ff.itp\" \n".format(GMX_FF_DIR))
    top_file.write("#include \"{}ff_b.itp\" \n".format(GMX_FF_DIR))
    top_file.write("; Include SPC water topology \n")
    top_file.write(";#include \"gromos53a6.ff/spc.itp\" \n") 
    top_file.write("#include \"{}spc_b.itp\" \n".format(GMX_FF_DIR))
    top_file.write("\n[ system ]\n")
    #top_file.write("Coarse-grained bilayer system\n")
    top_file.write("All-atom bilayer system\n")
    top_file.write("\n[ molecules ] \n") 
    

    return top_file

    # Using lipid system information, iterate through each molecule, printing out name and number of molecules 
    #for i, lipid_type in enumerate(lipid_system_info):
    #    n_molecule = int(lipid_type[1])
    #    molecule_name = lipid_type[0].name
    #    if n_molecule != 0:
    #        top_file.write("{:<10s}{:<10d}\n".format(molecule_name, n_molecule))

    ## Write out waters
    #top_file.write("{:<10s}{:<10d}\n".format(water().name, n_solvent))


def write_top_file_footer(top_file = None, n_solvent = 0):
    """ Generate topology file

    Parameters
    ----------
    top_file : File
        Topology file
    n_solvent : int
        Number of solvent beads in whole system

    Returns
    -------
    Updated topology file

    """

    # Just write out the number of solvents we have
    top_file.write("{:<10s}{:<10d}\n".format('SOL', n_solvent))
    return top_file

def solvate_bilayer(system = None, n_x = 8, n_y = 8, n_solvent_per_lipid = 5, water_spacing = 0.4, 
        res_index = 0, table_of_contents = None, lipid_atom_dict = None, atom_index = 0):
    """ Solvate the top and bottom parts of the bilayer, return water box

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
    lipid_atom_dict : OrderedDict()
        Dictionary whose values are mb.Compounds()s and values are a list of atom indices of that compound
    atom_index : int
        Counter for indexing atoms for lipid_atom_dict


    Returns
    -------
    solvated_system : mb.Compound()
        System with water solvating the outside of the bilayer
    water_box : mb.Box()
        Box object that accounts fot water molecules
    lipid_atom_dict : OrderedDict()
        Dictionary whose values are mb.Compounds()s and values are a list of atom indices of that compound
    atom_index : int
        Counter for indexing atoms for lipid_atom_dict

    """
    # Do some math to compute water spacing and box length/width/height
      # Construct 3D grid of water
    # Compute distances to translate such that water is either below or above bilayer
    # Add to table of contents file for post processing
    n_solvent_leaflet = n_x * n_y * n_solvent_per_lipid
    volume_solvent = 0.03 * n_solvent_leaflet # 0.03 nm3 is approximately the volume of one water molecule
    
    length = max(system.xyz[:,0])
    width = max(system.xyz[:,1])
    height = volume_solvent / (length * width)
    n_water_x = int(np.ceil(length/water_spacing))
    n_water_y = int(np.ceil(width/water_spacing))
    n_water_z = int(np.ceil(n_solvent_leaflet / (n_water_x * n_water_y)))
    cube = mb.Grid3DPattern(n_water_x, n_water_y, n_water_z)
    #cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    cube.scale( [ length, width, height])
    bot_water_list = cube.apply(H2O())
    bot_water_list = bot_water_list[ : n_solvent_leaflet]
    bot_water = mb.Compound()
    for compound in bot_water_list:
        bot_water.add(compound)
    highest_botwater = max(bot_water.xyz[:,2])
    lowest_botlipid = min(system.xyz[:,2])
    shift_botwater = abs(highest_botwater - lowest_botlipid) + 0.4
    mb.translate(bot_water, [0, 0, -1 * shift_botwater])
    # Add waters to table of contents
    #for i in range(n_x * n_y * n_solvent_per_lipid):
    for i in range(n_water_x * n_water_y * n_water_z):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, "HOH", H2O().n_particles))
        res_index += 1 
    
    # Construct 3D grid of water
    # Compute distances to translate such that water is either below or above bilayer
    # Add to table of contents file for post processing
    #cube = mb.Grid3DPattern(n_x, n_y, n_solvent_per_lipid)
    cube = mb.Grid3DPattern(n_water_x, n_water_y, n_water_z)
    #cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    cube.scale( [ length, width, height])
    top_water_list = cube.apply(H2O())
    top_water_list = top_water_list[ : n_solvent_leaflet]
    top_water = mb.Compound()
    for compound in top_water_list:
        top_water.add(compound)
    lowest_topwater = min(top_water.xyz[:,2])
    highest_toplipid = max(system.xyz[:,2])
    shift_topwater = abs(highest_toplipid - lowest_topwater) + 0.4
    mb.translate(top_water, [0, 0, shift_topwater])
    # Add waters to table of contents
    #for i in range(n_x * n_y * n_solvent_per_lipid):
    for i in range(n_water_x * n_water_y * n_water_z):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, "HOH", H2O().n_particles))
        res_index += 1

    system.add(bot_water)
    system.add(top_water)

    #waterbox = mb.Box(mins = [0,0,0], maxs = [max(top_water.xyz[:,0])+0.2, max(top_water.xyz[:,1])+0.2, max(top_water.xyz[:,2])+0.2])
    waterbox = mb.Box(mins = [0,0,0], maxs = [max(system.xyz[:,0])+0.2, max(system.xyz[:,1])+0.2, max(system.xyz[:,2])+0.2])

    
    # Add waters to lipid_atom_dict
    lipid_atom_dict['SOL'] = list(range(atom_index, atom_index + (2 * n_water_x * n_water_y * n_water_z * H2O().n_particles), 1))
    #atom_index += 2 * n_x * n_y * n_solvent_per_lipid * H2O().n_particles
    atom_index += 2 * n_water_x * n_water_y * n_water_z * H2O().n_particles
    n_solvent = 2 * n_water_x * n_water_y * n_water_z
    return system, waterbox, lipid_atom_dict, atom_index

def write_toc_file_box(table_of_contents = None, box = None):
    """ Write box information to table of contents file

    Parameters
    ----------
    table_of_contents : file
        Output file listing residue index, residue name, n_particles for tha residue

    box : mb.Box()
        System box


    """
    table_of_contents.write("{:<8.3f}{:<8.3f}{:<8.3f}\n".format(box.maxs[0], box.maxs[1], box.maxs[2]))

def write_ndx_file(filename = None, lipid_atom_dict = None):
    """ Write gromacs index file

    Parameters
    ---------
    filename : str
        filename to write .ndx file to
    lipid_atom_dict : OrderedDict()
        Dictionary whose values are mb.Compounds()s and values are a list of atom indices of that compound
    """
    ndx_file = open(filename+'.ndx','w')
    nonwater_string = ""
    water_string = ""
    system_string = ""
    
    for key in lipid_atom_dict.keys():
        # Join the indices into one big string
        #string_to_add = ' '.join([str(index) for index in lipid_atom_dict[key]])
        line = lipid_atom_dict[key]
        ndx_file.write(" [{}] \n".format(key))
        #chunk = []
        for i in range(0, len(line), 5):
            chunk = line[i:i+5]
            string_to_add = ' '.join([str(item) for item in chunk]) + '\n'


            # Print out that molecule group
            ndx_file.write("{}".format(string_to_add))

            # Get the water, nonwater, and system groups
             
            system_string += string_to_add
            if 'SOL' in key:
                water_string += string_to_add 
            else:
                nonwater_string += string_to_add
        
    ndx_file.write(" [ non-water ] \n")
    ndx_file.write(nonwater_string+"\n")
    ndx_file.write(" [ water ] \n")
    ndx_file.write(water_string+"\n")
    ndx_file.write(" [ System ] \n")
    ndx_file.write(system_string+"\n")
    ndx_file.close()

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "CG_bilayer", dest = "filename")
parser.add_option("-a", "--APL", action="store",type="float", default = 0.50, dest = "area_per_lipid")
parser.add_option("-r", "--rot", action="store", type ="float", default = 0.0, dest = "rotation")
parser.add_option("--DSPC", action="store",type="float", default = 1.0, dest = "DSPC_frac")
parser.add_option("--DPPC", action="store",type="float", default = 0.0, dest = "DPPC_frac")
parser.add_option("--acd16", action="store",type="float", default = 0.0, dest = "acd16_frac")
parser.add_option("--acd22", action="store",type="float", default = 0.0, dest = "acd22_frac")
parser.add_option("--alc12", action="store",type="float", default = 0.0,  dest = "alc12_frac")
parser.add_option("--alc14", action="store",type="float", default = 0.0, dest = "alc14_frac")
parser.add_option("--alc16", action="store",type="float", default = 0.0, dest = "alc16_frac")
parser.add_option("--alc18", action="store",type="float", default = 0.0, dest = "alc18_frac")
parser.add_option("--alc20", action="store",type="float", default = 0.0, dest = "alc20_frac")
parser.add_option("--alc22", action="store",type="float", default = 0.0, dest = "alc22_frac")
parser.add_option("--alc24", action="store",type="float", default = 0.0, dest = "alc24_frac")
parser.add_option("--acd24", action="store",type="float", default = 0.0, dest = "acd24_frac")
parser.add_option("--ISIS", action="store",type="float", default = 0.0, dest = "ISIS_frac")
parser.add_option("--SS", action="store",type="float", default = 0.0, dest = "SS_frac")
parser.add_option("--CHOL", action="store",type="float", default = 0.0, dest = "CHOL_frac")
parser.add_option("--PMEA", action="store",type="float", default = 0.0, dest = "PMEA_frac")
<<<<<<< HEAD
parser.add_option("--explicit", action="store_true" dest = "explicit")
=======
parser.add_option("--explicit", action="store_true", dest = "explicit")
>>>>>>> solvate_acd24
#parser.add_option("--water", action="store",type="float", default = 0.0, dest = "Water_frac")
(options, args) = parser.parse_args()


filename = options.filename
tilt_angle = options.rotation * np.pi/180
area_per_lipid = options.area_per_lipid
spacing = np.sqrt(area_per_lipid)

# Write out initial parameters
outfile = open((options.filename + 'initparam.txt'),'w')
outfile.write('Initial APL: {}\n'.format(options.area_per_lipid))
outfile.write('Initial Tilt: {}\n'.format(options.rotation))
outfile.close()


# Default parameters, less likely to be changed
n_x = 8#8
n_y = 8 #8
n_lipid = 2 * n_x * n_y
n_solvent_per_lipid = 20#5 # This is usually 20 waters per molecule, but a water bead is 4 waters
n_solvent = n_lipid * n_solvent_per_lipid
random_z_displacement = 0.3



# Todo: Obtain prototypes for more bilayer molecules
"""
lipid_system_info = [(DSPC(), np.ceil(n_lipid * options.DSPC_frac), 3.2),
                      (DPPC(), np.ceil(n_lipid * options.DPPC_frac), 2.6)
                      (acd16(), np.floor(n_lipid * options.acd16_frac), 4.4),
                      (acd22(), np.floor(n_lipid * options.acd22_frac), 3.2),
                      (alc12(), np.floor(n_lipid * options.Calc12_frac), 2.6),
                      (alc14(), np.floor(n_lipid * options.Calc14_frac), 2.8),
                      (alc16(), np.floor(n_lipid * options.Calc16_frac), 2.6),
                      (alc18(), np.floor(n_lipid * options.Calc18_frac), 2.4),
                      (alc20(), np.floor(n_lipid * options.Calc20_frac), 3.0),
                      (alc22(), np.floor(n_lipid * options.Calc22_frac), 3.2),
                      (alc24(), np.floor(n_lipid * options.Calc24_frac), 3.2),
                      (ISIS(), np.floor(n_lipid * options.ISIS_frac), 3.0),
                      (CHOL(), np.floor(n_lipid * options.CHOL_frac), 4.0)] 
                      """
"""
lipid_system_info = [(DSPC(), int(round(n_lipid * options.DSPC_frac)), 0.0),
                      (DPPC(), int(round(n_lipid * options.DPPC_frac)), 0.3),
                      (acd16(), int(round(n_lipid * options.acd16_frac)), 1.4),
                      (acd22(), int(round(n_lipid * options.acd22_frac)), 1.0),
                      (alc12(), int(round(n_lipid * options.Calc12_frac)), 1.4),
                      (alc14(), int(round(n_lipid * options.Calc14_frac)), 1.2),
                      (alc16(), int(round(n_lipid * options.Calc16_frac)), 1.0),
                      (alc18(), int(round(n_lipid * options.Calc18_frac)), 0.8),
                      (alc20(), int(round(n_lipid * options.Calc20_frac)), 0.6),
                      (alc22(), int(round(n_lipid * options.Calc22_frac)), 0.6),
                      (alc24(), int(round(n_lipid * options.Calc24_frac)), 0.6),
                      (ISIS(), int(round(n_lipid * options.ISIS_frac)), 0.4),
                      (CHOL(), int(round(n_lipid * options.CHOL_frac)), 0.8)] 
"""
# Hard coded lipid system info for complex systems
# This is for 50 WVTR
#lipid_system_info = [(DSPC(), 14, 0.0),
#                      (DPPC(), 0, -0.3),
#                      (acd16(), 0, -0.5),
#                      (acd22(), 22, -1.0),
#                      (alc12(), 0, -1.4),
#                      (alc14(), 0, -1.2),
#                      (alc16(), 48, -0.5),
#                      (alc18(), 0, -0.8),
#                      (alc20(), 0, -0.3),
#                      (alc22(), 30, -0.3),
#                      (alc24(), 0, -0.6),
#                      (ISIS(), 14, -2.5),
#                      (CHOL(), 0, -0.8)] 
#
#lipid_system_info = [(DSPC(), 14, 0.0),
#                      (DPPC(), 0, -0.3),
#                      (acd16(), 0, -0.5),
#                      (acd22(), 22, -1.0),
#                      (alc12(), 0, -1.4),
#                      (alc14(), 0, -1.2),
#                      (alc16(), 48, -0.5),
#                      (alc18(), 0, -0.8),
#                      (alc20(), 0, -0.3),
#                      (alc22(), 30, -0.3),
#                      (alc24(), 0, -0.6),
#                      (ISIS(), 14, -2.5),
#                      (CHOL(), 0, -0.8)] 

# For absolute numbers
if options.explicit:
    lipid_system_info = [(DSPC(), options.DSPC_frac, 0.0),
                       (DPPC(), options.DPPC_frac, -0.3),
                       (acd16(), options.acd16_frac, -1.4),
                       (acd22(), options.acd22_frac, -1.1),
                       (alc12(), options.alc12_frac, -1.4),
                       (alc14(), options.alc14_frac, -1.4),
                       (alc16(), options.alc16_frac, -1.4),
                       (alc18(), options.alc18_frac, -1.3),
                       (alc20(), options.alc20_frac, -1.2),
                       (alc22(), options.alc22_frac, -1.1),
                       (alc24(), options.alc24_frac, -1.0),
                       (acd24(), options.acd24_frac, -1.0),
                       (ISIS(), options.ISIS_frac, -3.0),
                       (CHOL(), options.CHOL_frac, -0.8)] 

# For doing fractions
else:
    lipid_system_info = [(DSPC(), np.ceil(n_lipid*options.DSPC_frac), 0.0),
                       (DPPC(), np.ceil(n_lipid*options.DPPC_frac), -0.3),
                       (acd16(), np.floor(n_lipid*options.acd16_frac), -1.4),
                       (acd22(), np.floor(n_lipid*options.acd22_frac), -1.1),
                       (alc12(), np.floor(n_lipid*options.alc12_frac), -1.4),
                       (alc14(), np.floor(n_lipid*options.alc14_frac), -1.4),
                       (alc16(), np.floor(n_lipid*options.alc16_frac), -1.4),
                       (alc18(), np.floor(n_lipid*options.alc18_frac), -1.3),
                       (alc20(), np.floor(n_lipid*options.alc20_frac), -1.2),
                       (alc22(), np.floor(n_lipid*options.alc22_frac), -1.1),
                       (alc24(), np.floor(n_lipid*options.alc24_frac), -1.0),
                       (acd24(), np.floor(n_lipid*options.acd24_frac), -1.0),
                       (ISIS(), np.floor(n_lipid*options.ISIS_frac), -3.0),
                       (CHOL(), np.floor(n_lipid*options.CHOL_frac), -0.8)] 

                      

#lipid_system_info = [(DSPC(), np.ceil(n_lipid * options.DSPC_frac), 0.0), #was 3.2
                     #(Calc12(), np.floor(n_lipid * options.Calc12_frac), 1.4)] #was 2.6
if(sum([lipid[1] for lipid in lipid_system_info]) != n_lipid):
    sys.exit("System setup error: number of components does not match layer size")

lipid_atom_dict = OrderedDict()
# Table of contents to externally store residue indices and names for post setup
table_of_contents = open(filename+'.dat', 'w')
res_index = 1
atom_index = 1

# Write topology file
print("Writing <{0}> ...".format(filename))
top_file = write_top_file_header(filename = filename, lipid_system_info = lipid_system_info)

# Generate bottom layer randomly
bot_layer, res_index, lipid_atom_dict, atom_index  = new_make_layer(n_x = n_x, n_y = n_y, lipid_system_info = lipid_system_info, 
        tilt_angle = tilt_angle, spacing = spacing, layer_shift = 3.2,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement, 
        top_file = top_file, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)

# Generate the top layer randomly
top_layer, res_index, lipid_atom_dict, atom_index  = new_make_layer(n_x = n_x, n_y = n_y, lipid_system_info = lipid_system_info, 
        tilt_angle = tilt_angle, spacing = spacing, layer_shift = 0,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement, 
        top_file = top_file, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)
# Rotate bottom layer to form bilayer
mb.spin_y(top_layer, theta=np.pi)

# Create system class that includes top and bottom layers
system = mb.Compound()
system.add(bot_layer)
system.add(top_layer)

# Solvate system, get new box
system, box, lipid_atom_dict, atom_index  = solvate_bilayer(system = system, n_x = n_x, n_y = n_y, n_solvent_per_lipid = n_solvent_per_lipid, 
        res_index = res_index, table_of_contents = table_of_contents, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)
top_file = write_top_file_footer(top_file = top_file, n_solvent = n_solvent)

"""
#new stuff
systembox=system._gen_box()
# Bilayer needs to be shifted to fit inside the box for gromacs (can't have any negative coordinates)
mb.translate(system, [np.abs(systembox.mins[0]), np.abs(systembox.mins[1]), np.abs(systembox.mins[2])])

# Shrink the x and y dimensions of the box to prevent packmol from shoving atoms in weird locations
cross_area = (n_x*spacing-1) * (n_y*spacing-1)
n_solvent_per_layer = 1280
leaflet_water_z = n_solvent_per_layer * 2.66602458e-2/cross_area
# Redefine the box to account for shifted coordinates and add extra space in z-direcction
smallbox = mb.Box(mins = [0, 0, systembox.mins[2] - leaflet_water_z], 
        maxs = [n_x * spacing, n_y * spacing, systembox.maxs[2] + leaflet_water_z])
# shift everything to be in box
min_z_shift = min(system.xyz[:,2])
mb.translate(system, [0.1, 0.1, -1 * min_z_shift])
smallbox.maxs[2] += -1*min_z_shift
# The solvate box will have a reduced cross-sectional area to prevent waters in strange locations
solvatebox = mb.Box(mins = [0.5, 0.5, smallbox.mins[2]],
    maxs = [smallbox.maxs[0] - 0.5, smallbox.maxs[1] - 0.5, smallbox.maxs[2]])


# Solvate system
system = mb.solvate(system, H2O(), 2*n_solvent_per_layer, solvatebox)

# Bilayer needs to be shifted to fit inside the box for gromacs (can't have any negative coordinates)
mb.translate(system, [np.abs(smallbox.mins[0]), np.abs(smallbox.mins[1]), np.abs(smallbox.mins[2])])

# Shift everything to positive z and off x and y axes
min_z_shift = min(system.xyz[:,2])
mb.translate(system, [0.1, 0.1, -1 * min_z_shift])
smallbox.maxs[2] += -1*min_z_shift

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    system.save(filename + '.gro', box =smallbox,overwrite=True)


# Add waters to lipid atom dict for index file
lipid_atom_dict['SOL'] = list(range(atom_index, atom_index + (2 * n_x * n_y * n_solvent_per_lipid * H2O().n_particles), 1))
atom_index += 2 * n_x * n_y * n_solvent_per_lipid * H2O().n_particles
write_ndx_file(filename = filename, lipid_atom_dict = lipid_atom_dict)


# Add waters to table of contents
for i in range(2* n_x * n_y * n_solvent_per_lipid):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, 'HOH', H2O().n_particles))
        res_index += 1
write_toc_file_box(table_of_contents = table_of_contents, box = smallbox)
# end new stuff
"""

# Script writer stuff
scriptWriter = scriptWriter("{}".format(options.filename)) 
scriptWriter.write_Titan_script(STrun = True)
scriptWriter.write_Titan_script(MDrun = True)
scriptWriter.write_Cori_script(STrun = True)
scriptWriter.write_Cori_script(MDrun = True)
scriptWriter.write_Rahman_script(STrun = True)
scriptWriter.write_Rahman_script(MDrun = True)


# Shift everything to positive z and off x and y axes
min_z_shift = min(system.xyz[:,2])
mb.translate(system, [0.1, 0.1, -1 * min_z_shift])
box.maxs[2] = max(system.xyz[:,2]) 

# Write to table of contents
write_toc_file_box(table_of_contents = table_of_contents, box = box)


# Write gro file suppressing mbuild warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    system.save(filename + '.gro', box =box,overwrite=True)
    #system.save(filename + '.gsd')
    #system.save(filename + '.xml')
table_of_contents.close()

# Write to an index file
write_ndx_file(filename = filename, lipid_atom_dict = lipid_atom_dict)

# Write job scripts
#thing1 = scriptWriter('{}'.format(filename))
#thing1.write_Rahman_script(MDrun = True)


