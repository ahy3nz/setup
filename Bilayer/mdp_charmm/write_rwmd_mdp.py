import os
import pdb
import random
import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.style.use('bmh')

def main():
    # First 25 ns is just the maximum range of heating
    cooling_rate = 1000 # 1 K every 1000 ps
    tc_groups = ['DSPC', 'alc12', 'water']
    t_pairs = [[305,575], [305,375], [305,505]]
    total_cooling_time = (t_pairs[0][1] - t_pairs[0][0])*cooling_rate 
    with open("heating_phase.mdp", 'w')as f:
        _write_body(f, ref_temp=305, n_steps=25000000, tc_groups=tc_groups)
        times, all_temps, t_pairs = _generate_heating_phase(f, t_pairs=t_pairs, 
                interval=5, time_anneal=25000, 
                dtemp=40, water_thermostat_style='plateau')

        _write_annealing_lines(f, times, all_temps)

    
    # Cooling phase is a series of 25ns simulations
    for i, t_start in enumerate(np.arange(25000, 25000 + total_cooling_time, 20000)):
        print("Writing cooling_phase{}.mdp".format(i))
        print("Temperature Pairs: {}".format(t_pairs))
        with open("cooling_phase{}.mdp".format(i),'w') as f:
            _write_body(f, ref_temp=305, n_steps=20000000, t_init=t_start, 
                    tc_groups=tc_groups)
            times, all_temps, t_pairs = _generate_cooling_phase(t_pairs=t_pairs,
                duration=20000, interval=5, cooling_rate=cooling_rate, 
                time_start=t_start, 
                water_thermostat_style='plateau')
            _write_annealing_lines(f, times, all_temps)



def _write_body(f, ref_temp=305, n_steps=50000000,t_init=0,tc_groups=None):
    tc_grps_string = " ".join([thing for thing in tc_groups])
    tau_t_string = " ".join(["0.4" for _ in range(len(tc_groups))])
    ref_t_string = " ".join(["{:5.0f}".format(ref_temp) for _ in range(len(tc_groups))])
    annealing_string = " ".join(["single" for _ in range(len(tc_groups))])
    f.write(""" title                       = RWMD
; Run parameters
integrator                  = md
nsteps                      = {n_steps}     ; 1fs/step * 5e7 steps = 5e7 fs = 50 ns
dt                          = 0.001
tinit                       = {t_init}

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
tc-grps                     = {tc_grps_string}    
tau_t                       = {tau_t_string}
ref_t                       = {ref_t_string}

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
annealing                   = {annealing_string}
""".format(**locals()))

def _generate_heating_phase(f,t_pairs=None,
        interval=5, time_anneal=25000, final_time=50000, dtemp=20,
        water_thermostat_style=None):
    """Generate annealing points
    f : File to write to 
    t_pairs : list of pairs,
        In each pair, (low temp, high temp). Convention
        is that the first ordered pair is the DSPC and the last 
        ordered pair is the water 
        """
    n_points = int(time_anneal/interval)

    #all_temps = [[] for _ in range(len(temp_pairs))]
    all_temps = np.zeros((len(t_pairs),n_points))

    # Initialize all temperatuers to their low temperature
    all_temps[:,0] = [temp[0] for temp in t_pairs]

    times = [0]

    # Generate scaling factors based on temperatre pairs
    if water_thermostat_style == 'proportional':
        scaling_factors = np.zeros(len(t_pairs))
        for i, t_pair in enumerate(t_pairs): 
            scaling_factors[i] = (t_pairs[i][1] - t_pairs[i][0]) / \
                (t_pairs[0][1] - t_pairs[0][0])
    

    possible_temps = np.arange(t_pairs[0][0], t_pairs[0][1]+1, dtemp)
    # Construct the histogram
    histogram = {val: 0 for val in possible_temps}
    histogram[all_temps[0,0]] += 1
    for i in np.arange(1,n_points):
        # new stuff here
        time = i*interval
        times.append(time)

        # Randomly pick a temperature
        new_temps = np.zeros(len(t_pairs))
        candidate_temp = random.choice(possible_temps)

        # Compare it to a criteria
        if histogram[all_temps[0,i-1]] >= histogram[candidate_temp]:
            new_temp = candidate_temp
            # Adjust temperature step for each tc group
            for j, tc_group in enumerate(t_pairs):
                if water_thermostat_style == 'proportional':
                    primary_temp_change = new_temp - t_pairs[0][0]
                    other_temp_change = scaling_factors[j] * primary_temp_change 
                    all_temps[j, i] = tc_group[0] + other_temp_change
            
                if water_thermostat_style == 'plateau':
                    if candidate_temp > tc_group[1]:
                        all_temps[j,i] = tc_group[1]
                    else:
                        all_temps[j,i] = candidate_temp
                    
        else:
            all_temps[:,i] = all_temps[:, i-1]

        histogram[all_temps[0,i]] += 1
        

    ### Saving stuff for reference
    with open('rwmd_histogram.dat', 'w') as thing:
        for k, v in histogram.items():
            thing.write("{}\t{}\n".format(k,v))
    #np.savetxt('timeseries.dat', np.column_stack((times, temps)))
    #### 

    #return times, temps, w_temps, (nonw_temp_low, nonw_temp_high), (w_temp_low, w_temp_high)
    return times, all_temps, t_pairs

    

def _generate_cooling_phase(time_start=25000, duration=50000, interval=5, 
        cooling_rate=500,
        t_pairs=None, 
        water_thermostat_style=None):
    """
    Change temperatures every 5 ps
    Reduce ceiling 2 K every ns
    Reduce ceiling 1 K every 500ps
    duration : int
        Length of cooling phase (ps)
    interval :
        duration of particular temperature before attempting jumpg
    cooling_rate : float
        Every `cooling_rate` ps, drop the RWMD temperature ceiling by 1 K
    water_thermostat_style : str
        Specifies how to thermstat water 
        'plateau' for water to plateau if thermostat exceeds a temperature
        'proportional' for water temp changes to be proportional to lipid temp

    """

    # This is a check to make sure temperatures aren't 
    # excessively higher than the primary temps
    n_points = int(np.floor(duration/interval))
    for i, t_pair in enumerate(t_pairs):
        if t_pair[1] > t_pairs[0][1]:
            t_pairs[i][1] = t_pairs[0][1]
        if t_pair[0] > t_pairs[0][0]:
            t_pairs[i][0] = t_pairs[0][0]
    
    if water_thermostat_style == 'proportional':
        scaling_factors = np.zeros(len(t_pairs))
        for i, t_pair in enumerate(t_pairs): 
            scaling_factors[i] = (t_pairs[i][1] - t_pairs[i][0]) / \
                (t_pairs[0][1] - t_pairs[0][0])

    times = [time_start]

    all_temps = np.zeros((len(t_pairs),n_points+1))

    # Initialize all temperatuers to their high temperature
    all_temps[:,0] = [temp[1] for temp in t_pairs]


    for i,time in enumerate(np.arange(time_start+interval, time_start+duration+interval, interval)):
        times.append(time)
        # If we've hit a cooling step, reduce the temperature range
        # Ensure that there is still a range of temperatures to choose from
        if time % cooling_rate == 0 and t_pairs[0][1] > t_pairs[0][0] + 1:
            t_pairs[0][1] -= 1

        new_temp = np.random.randint(t_pairs[0][0], t_pairs[0][1])

        for j, tc_group in enumerate(t_pairs):
            if water_thermostat_style == 'proportional':
                primary_temp_change = new_temp - t_pairs[0][0]
                other_temp_change = scaling_factors[j] * primary_temp_change 
                all_temps[j, i+1] = tc_group[0] + other_temp_change
        
            if water_thermostat_style == 'plateau':
                if new_temp > tc_group[1]:
                    all_temps[j,i+1] = tc_group[1]
                else:
                    all_temps[j,i+1] = new_temp 

        
    return times, all_temps, t_pairs 


def _write_annealing_lines(f, times, all_temps):
    if len(times) == np.shape(all_temps)[1]:
        temp_string = ""
        small_time_string = "".join(["{} ".format(thing) for thing in times])
        for i in range(all_temps.shape[0]):
            temp_string += "".join(["{:3.0f} ".format(thing) for thing in all_temps[i,:]])


        n_points = "".join(["{} ".format(len(times)) for _ in range(all_temps.shape[0])])
        time_string = " ".join([small_time_string for _ in range(all_temps.shape[0])])


        f.write("annealing-npoints\t={0}\n".format(n_points))
        f.write("annealing-time\t={0}\n".format(time_string))
        f.write("annealing-temp\t={0}\n".format(temp_string))
    else:
        sys.exit("Error, annealing lines do not have matching dimensions")





if __name__ == "__main__":
    main()

    
