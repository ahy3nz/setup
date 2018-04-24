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
#atom_types=['P4', 'P3', 'Nda', 'Na', 'C2', 'C1', 'Qa', 'Q0']
atom_types= ['W', 'E1','C2', 'C3', 'PCP', 'PCN']

# FF directory
FF_dir = '/home/yangah/Programs/setup/FF/CG/'
bond_parameters_file = os.path.join(FF_dir,'bond_parameters_msibi.dat')
angle_parameters_file = os.path.join(FF_dir, 'angle_parameters_msibi.dat')

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
    #for x,z in itertools.combinations_with_replacement(atom_types,2):
        for y in atom_types: 
            if '{}-{}-{}'.format(x,y,z) in angle_parameter_labels:
                angle_triplet = angle_parameters.loc['{}-{}-{}'.format(x,y,z)]
                harmonic_angles.angle_coeff.set(angle_triplet.name, k=angle_triplet.force_constant, t0=angle_triplet.x0)
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(z,y,x), k=angle_triplet.force_constant, t0=angle_triplet.x0)
                print("Setting angle {}-{}-{}".format(x,y,z))
                
            elif '{}-{}-{}'.format(z,y,x) in angle_parameter_labels:
                angle_triplet = angle_parameters.loc['{}-{}-{}'.format(z,y,x)]
                harmonic_angles.angle_coeff.set(angle_triplet.name, k=angle_triplet.force_constant, t0=angle_triplet.x0)
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=angle_triplet.force_constant, t0=angle_triplet.x0)
                print("Setting angle {}-{}-{}".format(z,y,x))
            else:
                print("WARNING: {}-{}-{} not found in angle parameters, setting to zero".format(x,y,z))
                harmonic_angles.angle_coeff.set('{}-{}-{}'.format(x,y,z), k=0, t0=0)


def set_bonds(system, constraints=True):
    """ Set all bond parameters
    Use a pandas dataframe to load in files

    Parameters
    ---------
    system : hoomd system object
    constraints : bool, default True
        If True, specifyc hoomd constraints for that bond


        """
    bond_parameters = pd.read_csv(bond_parameters_file, sep='\t',index_col=0)
    bond_harmonic = hoomd.md.bond.harmonic(name="mybond")
    for x,y in itertools.product(atom_types,repeat=2):
    #for x,y in itertools.combinations_with_replacement(atom_types,2):
        if '{}-{}'.format(x,y) in bond_parameters.index.values.tolist():
            bond_pair = bond_parameters.loc['{}-{}'.format(x,y)]
            bond_harmonic.bond_coeff.set(bond_pair.name, k=bond_pair.force_constant,
                        r0=bond_pair.x0)
            bond_harmonic.bond_coeff.set('{}-{}'.format(y,x), 
                    k=bond_pair.force_constant, r0=bond_pair.x0)
            print("Setting bond {}-{}".format(x,y))
            if constraints:
                _set_constraints(x, y, bond_pair.x0, system)
        elif '{}-{}'.format(y,x) in bond_parameters.index.values.tolist():
            bond_pair = bond_parameters.loc['{}-{}'.format(y,x)]
            bond_harmonic.bond_coeff.set(bond_pair.name, k=bond_pair.force_constant,
                        r0=bond_pair.x0)
            bond_harmonic.bond_coeff.set('{}-{}'.format(x,y), 
                    k=bond_pair.force_constant, r0=bond_pair.x0)
            print("Setting bond {}-{}".format(y,x))
            if constraints:
                _set_constraints(y, x, bond_pair.x0, system)
        else:
            print("WARNING: {}-{} not found in bond parameters, setting to zero".format(x,y))
            bond_harmonic.bond_coeff.set('{}-{}'.format(x,y), k=0, r0=0)
            bond_harmonic.bond_coeff.set('{}-{}'.format(y,x), k=0, r0=0)
    if constraints:
        constraints = hoomd.md.constrain.distance()
        constraint.set_params(rel_tol=0.01)
        return constraints

        
        
def set_pairs(table_dir="lambda0", nl=None, table=None, ignore=[]):
    """ load pair-wise interaction table for non-bonded interactions

    Parameters
    ---------
    table_dir = string
        Directory in which all the pair potnetials are located
    nl = hoomd neighborlist
    ignore = list of tuples
        [(i, j), (a, b),...] of pairs to ignore setting interactions

    Notes
    -----
    C1 is a 3:1 alkyl carbon mapping
    C2 is a 2:1 alkyl carbon mapping
    C2 pairs can jut use the C1 pairs



    """
    #tabulated_potential = os.path.join(FF_dir, table_dir)
    tabulated_potential = table_dir
    tablelength = 121
    if table is None:
        table = hoomd.md.pair.table(width = tablelength, nlist=nl)
    for atomtype_i, atomtype_j in itertools.combinations_with_replacement(atom_types,2):
        real_atomtype_i = str(atomtype_i).strip()
        #if 'C2' in real_atomtype_i:
        #    virtual_atomtype_i = 'C1'
        #else:
        #    virtual_atomtype_i = real_atomtype_i

        real_atomtype_j = str(atomtype_j).strip()
        #if 'C2' in real_atomtype_j:
        #    virtual_atomtype_j = 'C1'
        #else:
        #    virtual_atomtype_j = real_atomtype_j

        #if not set([real_atomtype_i, real_atomtype_j]) <= set(ignore):
        if True:
            table_file = '{}/{}-{}.txt'.format(tabulated_potential, real_atomtype_i, real_atomtype_j)
            table.set_from_file(real_atomtype_i, real_atomtype_j, filename=table_file)
            table.set_from_file(real_atomtype_j, real_atomtype_i, filename=table_file)
            print("Setting pair {}-{}".format(real_atomtype_i, real_atomtype_j))
        #elif set([atomtype_i, atomtype_j])<=set(['C1','C1']):
        #    table_file = '{}/{}-{}.txt'.format(tabulated_potential, atomtype_i, atomtype_j)
        #    table.set_from_file(atomtype_i, atomtype_j, table_file)
        #    table.set_from_file(atomtype_j, atomtype_i, table_file)


        return table

def set_exclusions(nl=None):
    """ Set 1-2 and 1-3 exclusions
        
    """
    excludelist=['1-3', '1-2']
    nl.reset_exclusions(exclusions = excludelist)


def _set_constraints(type_a, type_b, distance, system):
    """ Set bond constraints

    Parameters
    ---------
    type_a : str
        atomtype
    type_b : str
        atomtype
    distance : float
        Distance for constraint
    system : hoomd system object

    Returns
    -------
    constraint : hoomd.md.constrain.distance()
        """
    for bond in system.bonds:
        if (system.particles[bond.a].type == type_a and \
            system.particles[bond.b].type == type_b) or \
           (system.particles[bond.a].type == type_b and \
            system.particles[bond.b].type == type_a):
                system.constraints.add(index_a, index_b, distance)
