import mbuild as mb
from collections import OrderedDict
import warnings
import pdb
import numpy as np
import mdtraj as mdtraj
import sys
from itertools import product
from optparse import OptionParser
from Prototypes_CG import *
from scriptWriter import *

GMX_FF_DIR = "/raid6/homes/ahy3nz/Programs/setup/FF/CG/"
HOOMD_FF="/raid6/homes/ahy3nz/Programs/setup/FF/CG/myforcefield.xml"
# Mapping martini atomtpyes to atom names that are of that type
TYPE_TO_NAME_DICT = {'_Q0': ('_NC3', ''), '_Qa': ('_PO4', ''),
        '_Na': ('_GL1', '_GL2'), '_Nda': ('_LOH', '_COH', '_AOH', '_LgOH'),
        '_P2': ('_MOH', '_SOH', '_BOH'), '_P4': ('_W', '_COO'), '_BP4':  ('_WF', ''),
        '_C1': ('_C1', '_C2', '_C3', '_C4', '_C5', '_C1A', '_C2A', '_C3A', '_C4A', '_C5A', '_C1B', '_C2B', '_C3B', '_C4B', '_C5B') }

def convert_name_to_type(particle):
    """ Convert atom name to its corresponding atomtype

    Parameters
    ----------
    particle : mb.Compound()
    
    Returns
    -------

    """

    old_name = particle.name
    for atomtype in TYPE_TO_NAME_DICT:
        corresponding_atoms = TYPE_TO_NAME_DICT[atomtype]
        for single_atom in corresponding_atoms:
            #print(old_name+', ' + single_atom)
            if str(old_name) == str(single_atom):
                particle.name = atomtype
        #if str(old_name) in corresponding_atoms:
            #particle.name = atomtype
            else:
                pass

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

            # Rename the children of the molecule accordingly(not working)
            for child in molecule_to_add.children:
                convert_name_to_type(child)

            # Set charges (for the phospholipids)(not working)
            #if 'PC' in molecule_to_add.name:
            #    molecule_to_add.children[0].charge = 0.000122
            #    molecule_to_add.children[1].charge = -0.000122

            # When the pdb files are loaded, the coordinates are angstroms
            # Convert them to nm for all mbuild operations
            #molecule_to_add.xyz /= 10

            # Apply tilt angle
            molecule_to_add.spin(tilt_angle, [0,1,0])

            # Randomly orient them
            molecule_to_add.spin((np.pi/2)*np.random.rand() - (np.pi/4), [0,0,1])

            # Apply z_offset
            z_offset = lipid_type[2]

            # Apply APL and z_offset to identify the position for the molecule in the grid
            # Also generate random x and y offsets (+/- 0.1 nm)
            rand_x, rand_y = 0.2*np.random.rand(2) - 0.1
            position = [i * spacing + rand_x, j * spacing + rand_y, z_offset + layer_shift +
                        (-1 * np.random.random() * random_z_displacement)]
            molecule_to_add.translate(position)

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
    top_file.write("#include \"{}martini_ff.itp\" \n".format(GMX_FF_DIR))
    top_file.write(";#include \"{}martini_ff_b.itp\" \n".format(GMX_FF_DIR))
    top_file.write("\n[ system ]\n")
    top_file.write("Coarse-grained bilayer system\n")
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
    top_file.write("{:<10s}{:<10d}\n".format(w().name, n_solvent))
    return top_file

def solvate_bilayer(system = None, n_x = 8, n_y = 8, n_solvent_per_lipid = 5, water_spacing = 0.3, 
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
    length = max(system.xyz[:,0])
    width = max(system.xyz[:,1])
    n_water_leaflet = n_x * n_y * n_solvent_per_lipid # mol
    n_water_x = int(np.ceil(length/water_spacing)) # number of waters along x-direction
    n_water_y = int(np.ceil(width/water_spacing)) # number of waters along y-direction
    n_water_z = int(np.ceil(n_water_leaflet/(n_water_x * n_water_y))) # number of waters along z-direction
    height = n_water_z * water_spacing

    # Construct 3D grid of water
    # Compute distances to translate such that water is either below or above bilayer
    # Add to table of contents file for post processing
    #cube = mb.Grid3DPattern(n_x, n_y, n_solvent_per_lipid)
    cube = mb.Grid3DPattern(n_water_x, n_water_y, n_water_z)
    cube.scale([length, width, height])
    #cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    bot_water_list = cube.apply(w())
    bot_water_list = bot_water_list[ : n_water_leaflet]
    bot_water = mb.Compound()
    for compound in bot_water_list:
        for child in compound:
            convert_name_to_type(child)
            bot_water.add(compound)
    highest_botwater = max(bot_water.xyz[:,2])
    lowest_botlipid = min(system.xyz[:,2])
    shift_botwater = abs(highest_botwater - lowest_botlipid) + 0.3
    bot_water.translate([0, 0, -1 * shift_botwater])
    # Add waters to table of contents
    #for i in range(n_x * n_y * n_solvent_per_lipid):
    for i in range(n_water_x * n_water_y * n_water_z):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, w().name, w().n_particles))
        res_index += 1 
    
    # Construct 3D grid of water
    # Compute distances to translate such that water is either below or above bilayer
    # Add to table of contents file for post processing
    #cube = mb.Grid3DPattern(n_x, n_y, n_solvent_per_lipid)
    #cube.scale( [ water_spacing * n_x, water_spacing * n_y, water_spacing * n_solvent_per_lipid])
    cube = mb.Grid3DPattern(n_water_x, n_water_y, n_water_z)
    cube.scale([length, width, height])
    top_water_list = cube.apply(w())
    top_water_list = top_water_list[ : n_water_leaflet]
    top_water = mb.Compound()
    for compound in top_water_list:
        for child in compound:
            convert_name_to_type(child)
            top_water.add(compound)
    lowest_topwater = min(top_water.xyz[:,2])
    highest_toplipid = max(system.xyz[:,2])
    shift_topwater = abs(highest_toplipid - lowest_topwater) + 0.3
    top_water.translate([0, 0, shift_topwater])
    # Add waters to table of contents
    #for i in range(n_x * n_y * n_solvent_per_lipid):
    for i in range(n_water_x * n_water_y * n_water_z):
        table_of_contents.write("{:<10d}{:<10s}{:<10d}\n".format(res_index, w().name, w().n_particles))
        res_index += 1

    system.add(bot_water)
    system.add(top_water)

    waterbox = mb.Box(mins = [0,0,0], maxs = [max(top_water.xyz[:,0])+0.2, max(top_water.xyz[:,1])+0.2, max(top_water.xyz[:,2])+0.2])

    
    # Add waters to lipid_atom_dict
    #lipid_atom_dict['w'] = list(range(atom_index, atom_index + (2 * n_x * n_y * n_solvent_per_lipid * w().n_particles), 1))
    #atom_index += 2 * n_x * n_y * n_solvent_per_lipid * w().n_particles
    lipid_atom_dict['w'] = list(range(atom_index, atom_index + (2 * n_water_x * n_water_y * n_water_z * w().n_particles), 1))
    atom_index += 2 * n_water_x * n_water_y * n_water_z * w().n_particles
    n_solvent = 2 * n_water_x * n_water_y * n_water_z
    return system, waterbox, lipid_atom_dict, atom_index, n_solvent

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
            if 'w' in key:
                water_string += string_to_add 
            else:
                nonwater_string += string_to_add
        
    ndx_file.write(" [ Lipids ] \n")
    ndx_file.write(nonwater_string+"\n")
    ndx_file.write(" [ Water ] \n")
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
parser.add_option("--C16FFA", action="store",type="float", default = 0.0, dest = "C16FFA_frac")
parser.add_option("--C22FFA", action="store",type="float", default = 0.0, dest = "C22FFA_frac")
parser.add_option("--C12OH", action="store",type="float", default = 0.0,  dest = "C12OH_frac")
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
n_x = 8 #8
n_y = 8 #8
#n_x = 1
#n_y = 1
n_lipid = 2 * n_x * n_y
n_solvent_per_lipid = 20#5 # This is usually 20 waters per molecule, but a water bead is 4 waters
#n_solvent_per_lipid = 2
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

lipid_system_info = [(DSPC(), np.ceil(n_lipid * options.DSPC_frac), 0.0), #was 3.2
                     (OH12(), np.floor(n_lipid * options.C12OH_frac), 1.4)] #was 2.6

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
        tilt_angle = tilt_angle, spacing = spacing, layer_shift = 0,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement, 
        top_file = top_file, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)

# Generate the top layer randomly
top_layer, res_index, lipid_atom_dict, atom_index  = new_make_layer(n_x = n_x, n_y = n_y, lipid_system_info = lipid_system_info, 
        tilt_angle = tilt_angle, spacing = spacing, layer_shift = 4.0,
        res_index = res_index, table_of_contents = table_of_contents, random_z_displacement = random_z_displacement, 
        top_file = top_file, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)
       
# Rotate bottom layer to form bilayer
top_layer.spin(np.pi, [0,1,0])

# Create system class that includes top and bottom layers
system = mb.Compound()
system.add(bot_layer)
system.add(top_layer)

# Solvate system, get new box
system, box, lipid_atom_dict, atom_index, n_solvent = solvate_bilayer(system = system, n_x = n_x, n_y = n_y, n_solvent_per_lipid = n_solvent_per_lipid, water_spacing = 0.8,
        res_index = res_index, table_of_contents = table_of_contents, lipid_atom_dict = lipid_atom_dict, atom_index = atom_index)
top_file = write_top_file_footer(top_file = top_file, n_solvent = n_solvent)

# Shift everything to positive z and off x and y axes
min_z_shift = min(system.xyz[:,2])
system.translate( [0.1, 0.1, -1 * min_z_shift])
box = system.boundingbox
box.maxs[2] = max(system.xyz[:,2]) + 0.1
box.maxs[1] = max(system.xyz[:,1]) + 0.1
box.maxs[0] = max(system.xyz[:,0]) + 0.1



# Write to table of contents
write_toc_file_box(table_of_contents = table_of_contents, box = box)


# Write gro file suppressing mbuild warnings
#with warnings.catch_warnings():
    #warnings.simplefilter("ignore")
system.save(filename + '.gro', box =box,overwrite=True)
    
# In gromacs coordiantes, the smallest coordiante is just 0,0,0 
# We can just shift everything by the max/2
# Before saving hoomdxml, center everything around 0,0,0
# Also need to divide by 10 because of a unit conversion (pass as ref_distance)

system.translate([-box.maxs[0]/2, -box.maxs[1]/2, -box.maxs[2]/2])
box = system.boundingbox
box.lengths = box.lengths+2.0
# IF passing reference distances
system.save(filename + 'withref.hoomdxml', ref_energy = 0.239, ref_distance = 10, box=box, forcefield_files=HOOMD_FF, overwrite=True)
system.save(filename + 'withref.gsd', ref_energy = 0.239, ref_distance = 10,box=box, forcefield_files=HOOMD_FF,overwrite=True)
# For scaling by hand
system.xyz /= 10
box.lengths /=10
system.save(filename + 'scaled.hoomdxml', box=box, forcefield_files=HOOMD_FF, overwrite=True)
system.save(filename + 'scaled.gsd', box=box, forcefield_files=HOOMD_FF,overwrite=True)
table_of_contents.close()

# Write to an index file
write_ndx_file(filename = filename, lipid_atom_dict = lipid_atom_dict)

# Write job scripts
thing1 = scriptWriter('{}'.format(filename))
thing1.write_Rahman_script(MDrun = True)


