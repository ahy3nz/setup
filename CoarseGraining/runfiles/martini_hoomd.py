import hoomd
import hoomd.md
import numpy as np
import os
import sys
import pdb
import itertools

"""
HOOMD operates in arbitrary, self-consistent units. Martini uses gromacs units

Unit System
-----------
Energy = kj mol-1
length = nm
mass = amu = 1.6605e-27 kg
force = kj mol-1 nm-1
pressure = kj mol-1 nm-3 = 16.6054 bar
time = ps = 1e-12 s
charge = e = 1.60e2-19 C


Reduced Units (unitless)
------------
T_reduced = k * T / energy_unit
charge_reduced = q / ((4 * pi * eps_0 * distance_unit * energy_unit) ** 2)

Constants
---------
vacuum permittivity (eps_0) = 8.854e-12 C2 N-1 m-2
Boltzmann constant (k_B) = 1.381e-23 J K-1

"""
# These are conversion factors to convert from gromacs 
# units (kJ mol-1, nm, kj mol-1  nm-1) to table units
ENERGY_TABLE = 0.4184
DISTANCE_TABLE = 0.6
FORCE_TABLE = ENERGY_TABLE/DISTANCE_TABLE

# Martini atom types
atom_types = ['P5' , 'P4' , 'BP4', 'SP4', 'P3' , 'P2' ,  'P1' , 'SP1',
              'Nda', 'Nd' , 'SNd',  'Na' ,  'SNa', 'N0' ,
           'C5' ,  'SC5', 'C4' , 'SC4', 'C3' , 'SC3', 'C2' , 'AC2', 'SC2', 'C1' , 'C1a', 'AC1', 'SC1',
              'Qda', 'Qd', 'Qa' ,  'Q0' ]

# Martini interaction matrix (i.e. {}-{}: (epsilon, sigma))
interaction_atoms = ['Qda', 'Qd', 'Qa', 'Q0', 'P5', 'P4', 'P3', 'P2', 'P1', 'Nda', 'Nd', 'Na', 
    'N0', 'C5', 'C4', 'C3', 'C2', 'C1']
interaction_types = ["O", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"]
martini_epsilons = [5.6, 5, 4.5, 4, 3.5, 3.1, 2.7, 2.3, 2.0, 2.0]
martini_sigmas =   [0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.62]

# FF directory
FF_dir = '/raid6/homes/ahy3nz/Programs/setup/FF/CG/'
martini_interaction_matrix = os.path.join(FF_dir, "martini_interaction_matrix.dat")

def get_atom_types():
    """ Return a dictionary of all martini atom types"""
    return atom_types

def get_table_units():
    """ Return a tuple of table units"""
    return ENERGY_TABLE, DISTANCE_TABLE, FORCE_TABLE 


def set_harmonic_angles(timtails=False):
    """Set all angle parameters according to harmonic angles
    harmonic angle potentials were parametrized and fit to 
    cosine angle potentials provided by martini 

    Notes
    -----
    All angles are the same except for
    except for Qa-Na-Na, theta_0 = 120 deg 


    """
    # For using the force-fit harmonic parameters
    angle_kappa = 0.227 #kJ mol-1
    angle_kappa_unq = 18.58
    # For using the energy-fit harmonic parameters
    #angle_kappa = 0.135
    #angle_kappa_unq = 18.65

    # Tim parameters from Moore (2016) JPCB
    angle_kappa_timtail = 6.33 # This is in kcal mol-1 rad-2, 
    angle_kappa_timtail *= 4.184 # Convert to kJ mol-1 rad-2
    
    angle_theta0 = 180 #deg
    angle_theta0 = np.deg2rad(angle_theta0) # Convert to radians
    angle_theta0_unq = 120 #deg, for Qa-Na-Na
    angle_theta0_unq = np.deg2rad(angle_theta0_unq) # Convert to radians
    angle_theta0_timtail = 158
    angle_theta0_timtail = np.deg2rad(angle_theta0_timtail)
    harmonic_angles = hoomd.md.angle.harmonic()
    for x,z in itertools.product(atom_types,repeat=2):
        for y in atom_types: 
            # Angl-types containing 3xTim's C1 will be treated harmonically with his parameters
            #if (timtails) and ('C1' in (x,y,z)):
            if (timtails) and ('C1' in x) and ('C1' in y) and ('C1' in z) :
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k = angle_kappa_timtail, t0 = angle_theta0_timtail) 
                    
            # non-standard martini angle, Q0-Na-Na
            elif ('Na'== x and 'Na' == y and 'Qa' == z) or ('Qa' == x and 'Na' == y and 'Na' == z):
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k = angle_kappa_unq, t0 = angle_theta0_unq) 
            else:
                # Original angle
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k = angle_kappa, t0 = angle_theta0) 


def set_cosine_angles(timtails = False):
    """ Set hybrid angle potentials
    Tim's tails are treated with harmonic angles
    Martini tails are treated with cosine angles
    Set all angle parameters
     specify angles
     martini angles are V_angle(theta) = 0.5 * K_angle * (cos(theta) - cos(theta_0))**2
     all martini angles theta_0 = 180 deg and K_angle = 25 kJ*mol-1*nm-2
     except for Qa-Na-Na, theta_0 = 120 deg 


    Notes
    -----
    According to Martini FF, all angles are the same
    except for Qa-Na-Na, theta_0 = 120deg



    """
    # Martini parameters
    angle_kappa = 25 #kJ mol-1
    angle_theta0 = 180 #deg
    angle_theta0 = np.deg2rad(angle_theta0) # Convert to radians
    angle_theta0_unq = 120 #deg, for Qa-Na-Na
    angle_theta0_unq = np.deg2rad(angle_theta0_unq) # Convert to radians

    # Tim parameters from Moore (2016) JPCB
    angle_kappa_timtail = 6.33 # This is in kcal mol-1 rad-2, 
    angle_kappa_timtail *= 4.184 # Convert to kJ mol-1 rad-2
    angle_theta0_timtail = 158
    angle_theta0_timtail = np.deg2rad(angle_theta0_timtail)

    # Paramter v2, obtained from mapping AA sims to CG sims and getting this angle
    angle_kappa_2 = 25 #kJ mol-1
    angle_theta0_2 = 97.07
    angle_theta0_2 =  np.deg2rad(angle_theta0_2)

    cosine_angles = hoomd.md.angle.cosinesq()
    harmonic_angles = hoomd.md.angle.harmonic()

    # Iterate through and set all angle types for both the cosine and harmonic angle potentials
    for x,z in itertools.product(atom_types,repeat=2):
        for y in atom_types: 
            # Angletypes with C1 will be treated harmonically
            #if (timtails) and ('C1' in (x,y,z)):
            if (timtails) and ('C1' in x) and ('C1' in y ) and ('C1' in z):
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), 
                        k=angle_kappa_timtail, t0=angle_theta0_timtail) 
                cosine_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)
            # non-standard martini angle
            elif ('Na'== x and 'Na' == y and 'Qa' == z) or ('Qa' == x and 'Na' == y and 'Na' == z):
                cosine_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), 
                        k=angle_kappa, t0=angle_theta0_unq) 
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)
            # Angle v2, for C1a, Na, Na
            elif set(x,z) <= set('C1a', 'Na') and y == 'Na':
                cosine_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), 
                        k=angle_kappa_2, t0=angle_theta0_2)
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)

            # This angle is not actually parametrized, leave as 0
            elif set(x,z) <= set('C1', 'Na') and y =='Na':
                cosine_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)

            # Original martini angle
            else:
                cosine_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), 
                        k=angle_kappa, t0=angle_theta0)
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)


def set_bonds():
    """ set all bond parameters
    # specify bonds
    # martini bonds are V_bond(R) = 0.5 * K_bond * (R - R_bond)**2
    # all martini bonds have R_bond = 0.47 nm and K = 1250 kj*mol-1*nm-2
    # except for Na-Na bond, R_bond = 0.37nm
    """
    bond_k = 1250 #kJ mol-1 nm-2
    bond_r0 = 0.47 # nm
    bond_r0_unq = 0.37 # nm
    bond_harmonic = hoomd.md.bond.harmonic(name="mybond")
    for x,y in itertools.product(atom_types,repeat=2):
        # exception bond
        if 'Na' == x and 'Na' == y:
            bond_harmonic.bond_coeff.set('{}-{}'.format(x,y), k=bond_k, r0 = bond_r0_unq)
        # Standard bond
        else:
            bond_harmonic.bond_coeff.set('{}-{}'.format(x,y), k=bond_k, r0 = bond_r0)


def set_pairs(table_dir = "lambda0", nl = None):
    # load pair-wise interaction table for non-bonded interactions
    tabulated_potential = os.path.join(FF_dir, table_dir)
    tablelength = 121
    pair_table = hoomd.md.pair.table(width = tablelength, nlist=nl)
    for atomtype_i, atomtype_j in itertools.combinations_with_replacement(atom_types,2):
        atomtype_i = str(atomtype_i).strip()
        atomtype_j = str(atomtype_j).strip()
        table_file = '{}/{}-{}.txt'.format(tabulated_potential, atomtype_i, atomtype_j)
        pair_table.set_from_file(atomtype_i, atomtype_j, table_file)
    return pair_table

def set_exclusions(nl = None):
    """Tim's MSIBI tabulated potentials are built around 1,3 exclusions
    Martini only uses 1,2 exclusions
    Manually write ou the exclusiosn for the C1 tails
    So all C1-C1 interactions (the 1,3 C1-C1 interactions) need to be negelected
    """
    excludelist=['1-3', '1-2']
    nl.reset_exclusions(exclusions = excludelist)

def set_reactionfield_electrostatics(nl = None, r_cut = 1.1, eps_rf = 0.0, epsilon_r = 15):
    """ Martini uses reaction field electrostatics

    Parameters
    ---------
    nl : hoomd neighborlist
    r_cut : cutoff distance
    eps_rf : reaction field epsilon

    Notes
    -----
    Default r_cut, epsilon_r, and eps_rf pulled from martini parameters
    Set reactionfield pair-wise parameters according to Marrink 2007
    """
    # Load the interaction matrix file, excluding the first row/column
    #interaction_matrix_data = np.loadtxt(martini_interaction_matrix, skiprows = 1, usecols = np.arange(1,len(atom_types)))

    # Initialize the reaction field object
    reaction_field = hoomd.md.pair.reaction_field(r_cut = r_cut, nlist = nl)

    # Iterate and set the reaction field pair coefficients
    for x,y in itertools.combinations_with_replacement(interaction_atoms,2):
        # Turn the atomtype characters into an index for lists
        x_index = interaction_atoms.index(x)
        y_index = interaction_atoms.index(y)

        # Look at the interaction matrix to get the interaction classification
        #interaction_type = interaction_matrix_data[x_index,y_index]

        # Obtain the correct epsilon
        #epsilon = martini_epsilons[interaction_types.index(interaction_type)]

        # Set the reaction field parameters
        reaction_field.pair_coeff.set(x, y, epsilon = epsilon_r, eps_rf = eps_rf, use_charge = True)
