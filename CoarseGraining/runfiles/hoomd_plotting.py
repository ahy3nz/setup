import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

#filename = "log-output.log"
filename = "nvtnpt.log"

raw_data = np.loadtxt(filename, skiprows = 1)
cut_data = raw_data[1.5e6:,:]

#x_lim = [1.50e8, 1.80e8]
# timestep	potential_energy	temperature	lx	ly	lz
fig, axarray = plt.subplots(3,1 ,sharex=True)
axarray[0].plot(cut_data[:,0], cut_data[:,1])
axarray[0].set_title("Potential energy")
#axarray[0].set_ylim([-125000, -115000])
#axarray[0].set_xlim(x_lim)
axarray[1].plot(cut_data[:,0], cut_data[:,2])
axarray[1].set_title("Temperature")
#axarray[1].set_ylim([1.8, 2.7])
#axarray[1].set_xlim(x_lim)
axarray[2].plot(cut_data[:,0], cut_data[:,3]*cut_data[:,4]/64)
axarray[2].set_title("APL")
#axarray[2].set_ylim([0.5, 0.65])
#axarray[2].set_xlim(x_lim)
plt.savefig("CG_plots.png")
