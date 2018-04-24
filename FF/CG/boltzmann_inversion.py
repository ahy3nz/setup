import numpy as np
import scipy.interpolate
import glob
import pdb
from msibi.potentials import tail_correction, head_correction, alpha_array
from msibi.utils.smoothing import savitzky_golay


def _invert_rdf(rdf_i, kT):
    """ Invert RDF to get potnetial"""
    potential_i = -kT * np.log(rdf_i)
    return potential_i

def _linear_head(r, potential):
    """ Use linear extrapolation to fix
    potentials at small r"""
    last_bad_index = 0
    for i, val in enumerate(potential):
        if np.isnan(val) or np.isposinf(val) or np.isneginf(val):
            last_bad_index = i

    if potential[last_bad_index+1] < potential[last_bad_index+2]:
        f = scipy.interpolate.interp1d(r[last_bad_index+1 : last_bad_index +3], 
                potential[last_bad_index + 2 : last_bad_index+0:-1],
            fill_value='extrapolate')
    else:
        f = scipy.interpolate.interp1d(r[last_bad_index+1 : last_bad_index +3], 
            potential[last_bad_index + 1 : last_bad_index+3],
            fill_value='extrapolate')

    for i in np.arange(last_bad_index+1):
        potential[i] = f(r[i])
    return potential

rdfs = glob.glob('*-flhe_nvt.txt')
kT = 2.5
r_switch = 1.1
for rdf_file in rdfs:
    filename = "-".join(rdf_file.split('-')[0:2])
    print(filename)
    rdf = np.loadtxt(rdf_file) 
    dr = rdf[1,0] - rdf[0,0]
    potential = _invert_rdf(rdf[:,1], kT)
    potential = tail_correction(rdf[:,0], potential, r_switch)
    potential  = _linear_head(rdf[:,0], potential)
    forces = -1.0 * np.gradient(potential, dr)


    potential = savitzky_golay(potential,9,2,deriv=0,rate=1)
    np.savetxt(filename+".pot", np.column_stack((rdf[:,0], potential, forces)))


