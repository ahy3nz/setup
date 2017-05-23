import os
import sys
from scipy.optimize import curve_fit
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

"""
Convert the martini cosine angle potentials
to harmonic angle potnetials
using a fitting scheme for the forces
just outside the energy wells.
Consult the 2007 martini paper to
figure out the cosine angle parameters
"""
theta_0 = 180 #deg
k_cos = 25
theta_0 = np.deg2rad(theta_0)
theta_buffer = np.deg2rad(10)

def harmonic_force(theta, params):
    """ harmonic function to fit to
    
    Parameters
    ---------
    theta : independent variable, angle
        radians
    k_harm : harmonic angle constnat
    theta_0 : refernece angle

    Returns
    -------
    FORCE (not energy) evaluation 

    Notes:
    Hoomd angular energies are defined as
    V = 0.5 * k_harm * (theta - theta_0)**2
    Therefore,
    F = k_harm * (theta - theta_0)
    """
    k_harm = params
    force_eval = k_harm * (theta - theta_0)
    return force_eval

def harmonic_energy(theta, k_harm):
    """ Harmonic energy evaluation

    Parameters
    ---------
    theta : independent variable, angle
        radians
    k_harm : harmonic angle constnat
    theta_0 : refernece angle

    Returns
    -------
    Energy evaluation


    """
    return 0.5 * k_harm * ((theta - theta_0) ** 2)


# Define the independnet variable range, which is just 180 +/- 10 degrees
theta_range = np.linspace(theta_0 - theta_buffer, theta_0 + theta_buffer, num=100000)
# Define the dependent variables
cos_forces = k_cos * ( np.cos(theta_range) - np.cos(theta_0)) * (-np.sin(theta_range))
cos_energies =  0.5 * k_cos * ((np.cos(theta_range)- np.cos(theta_0))**2)

# DO the actual curve fitting
k_harm,covars = curve_fit(harmonic_force, theta_range, cos_forces) # Fit to forces
#k_harm,covars = curve_fit(harmonic_energy, theta_range, cos_energies) # Fit to energies

# generate some points
harmonic_forces = harmonic_force(theta_range, k_harm)
harmonic_energies = harmonic_energy(theta_range, k_harm)

fig, ax = plt.subplots(2,1)
ax[0].plot(np.rad2deg(theta_range), harmonic_forces, label = 'harmonic force')
ax[0].plot(np.rad2deg(theta_range), cos_forces, label='cosine force')
ax[0].legend()
ax[0].set_ylabel("Force")
ax[0].set_xlabel("Angle (degrees)")

ax[1].plot(np.rad2deg(theta_range), harmonic_energies, label = 'harmonic energy')
ax[1].plot(np.rad2deg(theta_range), cos_energies, label='cosine energy')
ax[1].legend()
ax[1].set_ylabel("Energy")
ax[1].set_xlabel("Angle (degrees)")
ax[0].set_title("k$^{harmonic}$="+str(k_harm))

plt.tight_layout()
plt.savefig('AngleFitting.png')
