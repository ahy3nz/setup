import os
import pdb
import random
import sys
import numpy as np

def _write_body(f, ref_temp=305):
    f.write(""" title                       = DSPC NPT Production Run
; Run parameters
integrator                  = md
nsteps                      = 50000000     ; 1fs/step * 5e7 steps = 5e7 fs = 50 ns
dt                          = 0.001

; Output control
nstxout                     = 0             ; Don't save coordinates 
nstvout                     = 0             ; Don't save velocities
nstenergy                   = 10000
nstlog                      = 10000
nstxtcout                   = 10000

;bond parameters
continuation                = yes
constraint_algorithm        = lincs
;constraints                 = all-bonds
constraints                 = none
lincs_iter                  = 1
lincs_order                 = 4

; Neighbor searching
cutoff-scheme               = Verlet
ns_type                     = grid
nstlist                     = 10
rcoulomb                    = 1.4
rvdw                        = 1.4

;Electrostatics
coulombtype                 = PME
pme_order                   = 4
fourierspacing              = 0.16

; Temperature coupling
tcoupl                      = nose-hoover
tc-grps                     = non-water water    
tau_t                       = 0.4   0.4     
ref_t                       = {ref_temp}   {ref_temp}

;Pressure coupling
pcoupl                      = Parrinello-Rahman
pcoupltype                  = semiisotropic
tau_p                       = 2.0           ; ps
ref_p                       = 1.0 1.0          ; bar   
compressibility             = 4.5e-5 4.5e-5
refcoord_scaling            = com

;PBC
pbc                         = xyz

;Dispersion correction
DispCorr                    =EnerPres

;Velocity generation
gen_vel                     = no

;Simulated annealing
annealing                   = single single
""".format(**locals()))

def _generate_annealing_points(f, temp_low=305, temp_high=560, interval=5, 
        time_anneal=25000, final_time=50000, water_temp=500, dtemp=20):
    """ Generate annealing points
    temp_low : lowest temperature, float
    temp_high : highest temperature, float
    interval : duration for each annealing interval, int (ps)
    time_anneal: duration of annealing, int (ps) 
    final_time : last time point such that the RWMD gradually
        brings temperature down to temp_low at final_time, int (ps)
    water_temp : temperature to keep water at, K 
    dtemp : temperature intervals, (K)
    """
    n_points = int(time_anneal/interval)
    times = [0]
    temps = [temp_low]
    last_temp = temp_low
    possible_temps = np.arange(temp_low, temp_high+1, dtemp)
    # Construct the histogram
    histogram = {val: 0 for val in possible_temps}
    histogram[last_temp] += 1
    for i in np.arange(1,n_points):
        # new stuff here
        time = i*interval
        times.append(time)

        # Randomly pick a temperature
        candidate_temp = random.choice(possible_temps)

        # Compare it to a criteria
        if histogram[last_temp] >= histogram[candidate_temp]:
            new_temp = candidate_temp
        else:
            new_temp = last_temp

        histogram[new_temp] += 1
        temps.append(new_temp)
        last_temp = new_temp

        ## ---old suf---
        ## If we're at the lowest temperature, only go up
        #if abs(last_temp - temp_low) < 0.1:
        #    temp_change = dtemp * np.random.randint(0,2)
        ## If we're at the highest temperature, only go down
        #elif abs(last_temp - temp_high) < 0.1:
        #    temp_change = dtemp * np.random.randint(-1,1)
        #else:
        #    temp_change = dtemp * np.random.randint(-1,2)
        #time = i*interval
        #new_temp = last_temp + temp_change
        #times.append(time)
        #temps.append(new_temp)
        #last_temp = new_temp
        ## -----end---

    ### Saving stuff for reference
    #with open('dictionary.dat', 'w') as thing:
    #    for k, v in histogram.items():
    #        thing.write("{}\t{}\n".format(k,v))
    #np.savetxt('timeseries.dat', np.column_stack((times, temps)))
    #### 

    temp_string = " ".join([str(temp) for temp in temps])
    temp_string += " " + str(temp_low)

    time_string = " ".join([str(time) for time in times])
    time_string += " " + str(final_time)
    f.write("annealing-npoints\t={0} 2 \n".format(n_points+1))
    f.write("annealing-time\t={0} 0 {1}\n".format(time_string, str(final_time)))
    f.write("annealing-temp\t={0} {1} {1}\n".format(temp_string, str(water_temp)))

if __name__ == "__main__":
    with open('RW.mdp','w') as f:
        _write_body(f, ref_temp=305)
        _generate_annealing_points(f, temp_low=305, temp_high=585, interval=5,
                time_anneal=25000, final_time=50000, water_temp=305, dtemp=28)
