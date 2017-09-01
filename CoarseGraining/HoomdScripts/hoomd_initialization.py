import hoomd
import hoomd.md
import numpy as np
import os
import pdb
import pandas as pd
import itertools

"""
HOOMD operates in arbitrary, self-consistent units. 
These simulations use gromacs units

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
# units (kJ mol-1, nm, kj mol-1  nm-1) to Tim table units
ENERGY_TABLE = 0.4184
DISTANCE_TABLE = 0.6
FORCE_TABLE = ENERGY_TABLE/DISTANCE_TABLE

# Atom types
atom_types=['P4', 'P3', 'Nda', 'Na', 'C1', 'Qa', 'Q0']

# FF directory
FF_dir = '/raid6/homes/ahy3nz/Programs/setup/FF/CG/'
bond_parameters_file = os.path.join(FF_dir,'bond_parameters.dat')
angle_parameters_file = os.path.join(FF_dir, 'angle_parameters.dat')

def get_atom_types():
    """ Return a dictionary of all atom types"""
    return atom_types

def get_table_units():
    """ Return a tuple of table units"""
    return ENERGY_TABLE, DISTANCE_TABLE, FORCE_TABLE 


def set_harmonic_angles():
    """ Set all harmonic angles 
    Use a pandas dataframe to load in files 

    """

    angle_parameters = pd.read_csv(angle_parameters_file, sep='\t',index_col=0)
    angle_parameter_labels = angle_parameters.index.values.tolist()
    harmonic_angles = hoomd.md.angle.harmonic()
    for x,z in itertools.product(atom_types,repeat=2):
        for y in atom_types: 
            if '{}-{}-{}'.format(x,y,z) in angle_parameter_labels:
                angle_triplet = angle_parameters.loc('{}-{}-{}'.format(x,y,z))
                harmonic_angles.angle_coeff.set(angle_triplet.name, k=angle_triplet.force_constant, t0=angle_triplet.x0)
                
            elif '{}-{}-{}'.format(z,y,x) in angle_parameter_labels:
                angle_triplet = angle_parameters.loc('{}-{}-{}'.format(z,y,x))
                harmonic_angles.angle_coeff.set(angle_triplet.name, k=angle_triplet.force_constant, t0=angle_triplet.x0)
            else:
                print("WARNING: {}-{}-{} not found in angle parameters, setting to zero".format(x,y))
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)


def set_bonds():
    """ Set all bond parameters
    Use a pandas dataframe to load in files


        """
    bond_parameters = pd.read_csv(bond_parameters_file, sep='\t',index_col=0)
    bond_harmonic = hoomd.md.bond.harmonic(name="mybond")
    for x,y in itertools.product(atom_types,repeat=2):
        if '{}-{}'.format(x,y) in bonding_parameters.index.values.tolist()
        try:
            bond_pair = bond_parameters.loc('{}-{}'.format(x,y))
            bond_harmonic.bond_coeff.set(bond_pair.name, k=bond_pair.force_constant,
                    r0=bond_pair.x0)
        except KeyError:
            print("WARNING: {}-{} not found in bond parameters, setting to zero".format(x,y))
            bond_harmonic.bond_coeff.set('{}-{}'.format(x,y), k=0, r0=0)
        
        
def set_pairs(table_dir = "lambda0", nl = None, table=None, ignore=None):
    # load pair-wise interaction table for non-bonded interactions
    tabulated_potential = os.path.join(FF_dir, table_dir)
    tablelength = 121
    if table is None:
        table = hoomd.md.pair.table(width = tablelength, nlist=nl)
    for atomtype_i, atomtype_j in itertools.combinations_with_replacement(atom_types,2):
        atomtype_i = str(atomtype_i).strip()
        atomtype_j = str(atomtype_j).strip()
        if not set([atomtype_i, atomtype_j]) <= set(ignore):
            table_file = '{}/{}-{}.txt'.format(tabulated_potential, atomtype_i, atomtype_j)
            table.set_from_file(atomtype_i, atomtype_j, table_file)
            table.set_from_file(atomtype_j, atomtype_i, table_file)
        elif set([atomtype_i, atomtype_j])<=set(['C1','C1']):
            table_file = '{}/{}-{}.txt'.format(tabulated_potential, atomtype_i, atomtype_j)
            table.set_from_file(atomtype_i, atomtype_j, table_file)
            table.set_from_file(atomtype_j, atomtype_i, table_file)

    return table

def set_exclusions(nl=None):
    """ Set 1-2 and 1-3 exclusions
        
    """
    excludelist=['1-3', '1-2']
    nl.reset_exclusions(exclusions = excludelist)


