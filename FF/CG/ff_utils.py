import os
import pdb
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from optparse import OptionParser
"""
Convert Martini (or any gromacs) LJ parameters 
to tabulated form for HOOMD simulation.
table_energy [=] 0.1 kcal/mol
table_dist [=] 6 angstrom
table_force [=] table_energy/table_dist
specify f, c6, c12
Martini cutoff at 1.1nm
"""
# These are conversion factors to convert from gromacs 
# units (kJ mol-1, nm, kj mol-1  nm-1) to table units
ENERGY_GMX_TABLE = 0.4184
DISTANCE_GMX_TABLE = 0.6
FORCE_GMX_TABLE = ENERGY_GMX_TABLE/DISTANCE_GMX_TABLE

def convert_unit(value=0, tag=None):
    """ Convert a gromacs unit to tabular unit

    Parameters
    ---------
    value : float
        value to be converted
    tag : str
        specify fundamental unit

    Returns
    -------
    converted_unit : float
        Value, converted to tabular units

"""
    if tag == "energy":
        converted_unit = value / ENERGY_GMX_TABLE
    elif tag == "force":
        converted_unit = value / FORCE_GMX_TABLE
    elif tag == "distance":
        converted_unit = value / DISTANCE_GMX_TABLE
    else:
        sys.exit("Specify a tag (energy, force, distance) for unit conversion")
    return converted_unit

def calc_LJ_energy(r = 1, c12 = 0, c6 = 0, rcut = 9999999):
    """ Compute energy according to LJ potential

    Parameters
    ----------
    r : float
        distance (nm)
    c12 : float
        c12 parameter (see section 4.1.1 in gmx manual)
    c6 : float
        c6 parameter (see section 4.1.1 in gmx manual)

    Returns
    -------
    energy : float
        energy (see section 4.1.1 in gmx manual)

    Notes
    -----
    This LJ energy calculation is based off gromacs units.
    See section 4.1.1 in the gromacs manual
    
    """
    # Exception for r = 0
    if r == 0 or r >= rcut:
        return 0
    r12 = r ** 12
    r6 = r ** 6
    energy = (c12/r12) - (c6/r6)
    return energy

def calc_LJ_force(r = 1, c12 = 0, c6 = 0):
    """ Compute force according to derivative of LJ potential

    Parameters
    ----------
    r : float
        distance (nm)
    c12 : float
        c12 parameter (see section 4.1.1 in gmx manual)
    c6 : float
        c6 parameter (see section 4.1.1 in gmx manual)

    Returns
    -------
    force : float
        force (see section 4.1.1 in gmx manual)

    Notes
    -----
    This LJ force calculation is based off gromacs units.
    See section 4.1.1 in the gromacs manual
    """
    # Exception for r = 0
    if r == 0:
        return 0
    r13 = r ** 13
    r7 = r ** 7
    force = ((12*c12/r13) - (6*c6/r7))
    return force


def calc_morse_energy(r, D, alpha, r0):
    """Morse pair potential. """
    return D * (np.exp(-2 * alpha * (r - r0)) - 2 * np.exp(-alpha * (r - r0)))

def calc_morse_force(r, D, alpha, r0):
    """ Morse pair force"""
    return D*(-2*alpha*np.exp(-alpha*(r-r0)) + 2*alpha*np.exp(-2*alpha*(r-r0)))

def LJ_to_morse(start_fit=0.3, end_fit=0.6, sigma=1, eps=1, output_plot=False):
    """ Take LJ parameters and fit to a morse potential
    
    Parameters
    ----------
    start_fit : float
        Fit to a particular region of the LJ potential
    end_fit : float
        Fit to a particular region of hte LJ potential

    Returns
    -------
    morse_params : dict

    """
    distance_vals = np.arange(start_fit, end_fit, 0.01)
    c6 = 4*eps*(sigma**6)
    c12 = 4*eps*(sigma**12)
    LJ_energies = [calc_LJ_energy(r=r, c6=c6, c12=c12) for r in distance_vals]

    params, covar = curve_fit(calc_morse_energy, distance_vals, LJ_energies, 
            bounds=[(0,0,0), (np.inf, np.inf, np.inf)])
    morse_params = {'D': params[0], 'alpha': params[1], 'r0': params[2]}


    return morse_params




def interpolate_r0(values = None, distances = None):
    """ Compute values at r=0 by interpolating
    
    Parameters
    ---------
    values = list
        List whose indices corresopnd to different distances
    distances = list
        List of distances (correspond to elements in values)

    Returns
    -------
    modified_values = list
        List whose indices corresopnd to different distances
"""
    modified_values = values
    dy = values[2] - values[1]
    dx = distances[2] - distances[1]
    slope = dy/dx
    dx_0 = distances[1] - distances[0]
    modified_values[0] = values[1] - (slope*dx_0)
    return modified_values

def gmx_to_sigeps(c6=0, c12=0, type_A=None, type_B=None):
    """ Compute sigma and epsilon from gromacs c6 and c12

    Parameters
    ----------
    c12 : float
        c12 parameter (see section 4.1.1 in gmx manual)
    c6 : float
        c6 parameter (see section 4.1.1 in gmx manual)
    type_A : str
        atom type
    type_B : str
        atom type

    Returns
    ------
    sig : float (nm)
        sigma, radius
    eps : float (kJ/mol)
        epsilon, depth of well
    """
    sig = (c12/c6)**(1/6)
    eps_1 = c12/(4*(sig**12))
    eps_2 = c6/(4*(sig**6))
    if round(eps_1, 4) != round(eps_2, 4):
        print("Numerical error with epsilon checking {}-{}".format(type_A, type_B))
    return sig, eps_1


def generate_morse_Table(output=None, plot=False, D=1, alpha=1, r0=1):
    """ Generated rabulated tabulated morse potentials

    """

    distances = np.linspace(0, 1.2, num = 121) 
    rcut = 1.1 
    energies = [calc_morse_energy(r, D, alpha, r0) for r in distances] 
    forces = [calc_morse_force(r, D, alpha, r0) for r in distances] 

    # Set energies to zero beyond cutoff
    for i, distance in enumerate(distances):
        if distance >= rcut:
            energies[i] = 0

        
    # Save to text file
    zipped =["{:.18e} {:.18e} {:.18e}".format(dist, ener, force) for dist,ener,force in zip(distances, energies, forces)] 
    np.savetxt('{}.txt'.format(output), zipped, fmt='%s')

    
    # Plot tabulated potentials
    if plot:
        f, axarray= plt.subplots(1,1)
        axarray.plot(distances, energies, label="energy")
        axarray.plot(distances, forces, label = "force")
        axarray.plot(distances, np.zeros(len(distances)))
        axarray.set_ylim([-10, 10])
        axarray.set_title("Gromacs")
        axarray.legend()
        
        plt.savefig(output+".png")
        plt.close()

def generate_Table(c6=0, c12=0, output=None, plot=False):
    """ Generated tabulated pair potential

    Parameters
    ---------
    c6 : float
        c6 term for gromacs
    c12 : float
        c12 term for gromacs
    output : string
        output filename
    plot : boolean
        Plot forces and energies
        """
    distances = np.linspace(0, 1.2, num = 121) 
    rcut = 1.1 
    energies = [calc_LJ_energy(r=r, c6=c6, c12=c12, rcut=rcut) for r in distances] 
    forces = [calc_LJ_force(r=r, c6=c6, c12=c12) for r in distances] 
    energies = interpolate_r0(values=energies, distances=distances) 
    forces = interpolate_r0(values=forces, distances=distances) 
        
    # Save to text file
    zipped =["{:.18e} {:.18e} {:.18e}".format(dist, ener, force) for dist,ener,force in zip(distances, energies, forces)] 
    np.savetxt('{}.txt'.format(output), zipped, fmt='%s')

    
    # Plot tabulated potentials
    if plot:
        f, axarray= plt.subplots(1,1)
        axarray.plot(distances, energies, label="energy")
        axarray.plot(distances, forces, label = "force")
        axarray.plot(distances, np.zeros(len(distances)))
        axarray.set_ylim([-10, 10])
        axarray.set_title("Gromacs")
        axarray.legend()
        
        plt.savefig(output+".png")
        plt.close()

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", action = "store", type = "string", dest = "filename")
    parser.add_option("--c6", action = "store", type = "float", default = 0.00, dest = "c6")
    parser.add_option("--c12", action = "store", type = "float", default = 0.00, dest = "c12")
    (options, args) = parser.parse_args()
    c6 = options.c6
    c12 = options.c12

    if not options.filename:
        sys.exit("specify filename")
    if options.c6 == 0.00:
        sys.exit("specify c6")
    if options.c12 == 0.00:
        sys.exit("specify c12")


    # Perform all calculations using gromacs units, convert at the end

