import glob
import os
import shutil
from scipy.optimize import curve_fit
import numpy as np
import itertools 

def scale_region(filename, zone, scale_factor, switch):
    potential = np.loadtxt(filename)
    step = potential[1,0] - potential[0,0]
    for i in range(potential.shape[0]):
        # Do full scaling if within the zone
        if zone[0] <= potential[i,0] <= zone[1] and potential[i,1] < 0:
            potential[i,1] *= scale_factor
        # Do partial scaling if outside zone, within switching distance
        elif ((min(abs(zone[0] - potential[i,0]), 
                abs(zone[1] - potential[i,0])) < switch) and 
                potential[i,1] < 0):
            dist = min(abs(zone[0] - potential[i,0]), abs(zone[1] - potential[i,0]))
            partial_scale = (1 - (dist/switch)) * scale_factor
            potential[i,1] *= partial_scale
        else:
            pass
    potential[:,2] = -1* np.gradient(potential[:,1], step)
    np.savetxt(filename, potential)

def line_func(x, m, b):
        return m*x + b

def shift_routine(file_string, shift_outward, dr=0.01):
    potential = np.loadtxt(file_string)
    new_r_vals = np.arange(0, shift_outward, step=dr)
    popt, pcov = curve_fit(line_func, potential[:5,0]+shift_outward, potential[:5,1])
    new_v_vals = line_func(new_r_vals, *popt)
    combined_v_vals = np.concatenate((new_v_vals, potential[:,1]))
    
    new_arr = np.column_stack((potential[:,0], combined_v_vals[0:121]))

    new_grad = -1 * np.gradient(new_arr[:,1], dr)
    np.savetxt(file_string , np.column_stack((new_arr, new_grad))[0:121])

def shift_inward(file_string, shift_inward, dr=0.01):
    potential = np.loadtxt(file_string)
    n_rows = int(shift_inward/dr)
    insert_lines = np.zeros((n_rows, ))
    combined_v_vals = np.concatenate(( potential[:,1], insert_lines))
    new_grad = -1 * np.gradient(combined_v_vals, dr)
    np.savetxt(file_string, np.column_stack((potential[:,0], 
        combined_v_vals[n_rows:], new_grad[n_rows:])))

def scale_potential(file_string, scaling_factor, dr=0.01):
    potential = np.loadtxt(file_string)
    new_V = potential[:,1] * scaling_factor
    new_grad = -1 * np.gradient(new_V, dr)
    np.savetxt(file_string, np.column_stack((potential[:,0],
        new_V, new_grad)))
