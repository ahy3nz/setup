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
    nonw_temp_high = 705
    nonw_temp_low = 305
    total_cooling_time = (nonw_temp_high - nonw_temp_low)*cooling_rate 
    with open("heating_phase.mdp", 'w')as f:
        _write_body(f, ref_temp=305, n_steps=25000000)
        times, temps, w_temps, nonw_temp_range, w_temp_range = _generate_heating_phase(f, nonw_temp_low=nonw_temp_low, nonw_temp_high=nonw_temp_high,
                interval=5, time_anneal=25000, w_temp_low=305, w_temp_high=505,
                dtemp=40, water_thermostat_style='plateau')

        _write_annealing_lines(f, times, temps, w_temps)

    
    # Cooling phase is a series of 25ns simulations
    for i, t_start in enumerate(np.arange(25000, 25000 + total_cooling_time, 25000)):
        print("Writing cooling_phase{}.mdp".format(i))
        print("NonW: {}-{}".format(nonw_temp_range[0], nonw_temp_range[1]))
        with open("cooling_phase{}.mdp".format(i),'w') as f:
            _write_body(f, ref_temp=305, n_steps=25000000, t_init=t_start)
            times, temps, w_temps, nonw_temp_range, w_temp_range = _generate_cooling_phase(nonw_temp_low=nonw_temp_range[0], nonw_temp_high=nonw_temp_range[1],
                duration=25000, interval=5, cooling_rate=cooling_rate, 
                time_start=t_start, 
                w_temp_high=w_temp_range[1], w_temp_low=w_temp_range[0],
                water_thermostat_style='plateau')
            _write_annealing_lines(f, times, temps, w_temps)



def _write_body(f, ref_temp=305, n_steps=50000000,t_init=0):
    f.write(""" title                       = DSPC NPT Production Run
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

def _generate_heating_phase(f, nonw_temp_low=305, nonw_temp_high=560, interval=5, 
        time_anneal=25000, final_time=50000, w_temp_low=500, w_temp_high=500,
        dtemp=20,
        water_thermostat_style=None):
    """ Generate annealing points
    nonw_temp_low : lowest temperature, float
    nonw_temp_high : highest temperature, float
    interval : duration for each annealing interval, int (ps)
    time_anneal: duration of annealing, int (ps) 
    final_time : last time point such that the RWMD gradually
        brings temperature down to nonw_temp_low at final_time, int (ps)
    w_temp_low : temperature to keep water at, K 
    w_temp_high : temperature to keep water at, K 
    dtemp : temperature intervals, (K)
    water_thermostat_style : str
        Specifies how to thermstat water 
        'plateau' for water to plateau if thermostat exceeds a temperature
        'proportional' for water temp changes to be proportional to lipid temp
        'constant' for water to maintain a single temperature
            """
    n_points = int(time_anneal/interval)
    times = [0]
    temps = [nonw_temp_low]
    if water_thermostat_style == 'proportional':
        w_temps= [w_temp_low]
        last_w_temp = w_temp_low
        scaling_factor = (w_temp_high - w_temp_low)/(nonw_temp_high - nonw_temp_low)
    if water_thermostat_style == 'plateau':
        w_temps = [w_temp_low]
        last_w_temp = w_temp_low


    last_temp = nonw_temp_low
    possible_temps = np.arange(nonw_temp_low, nonw_temp_high+1, dtemp)
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
            if water_thermostat_style == 'proportional':
                nonw_temp_change = new_temp - nonw_temp_low
                w_temp_change = scaling_factor * nonw_temp_change 
                new_w_temp = w_temp_low + w_temp_change
            if water_thermostat_style == 'plateau':
                if candidate_temp > w_temp_high:
                    new_w_temp = w_temp_high
                else:
                    new_w_temp = candidate_temp

        else:
            new_temp = last_temp
            if water_thermostat_style == 'proportional':
                new_w_temp = last_w_temp
            if water_thermostat_style == 'plateau':
                new_w_temp = last_w_temp

        histogram[new_temp] += 1
        temps.append(new_temp)
        if water_thermostat_style == 'proportional':
            w_temps.append(new_w_temp)
            last_w_temp = new_w_temp
        if water_thermostat_style == 'plateau':
            w_temps.append(new_w_temp)
            last_w_temp = new_w_temp

        last_temp = new_temp


    ### Saving stuff for reference
    #with open('dictionary.dat', 'w') as thing:
    #    for k, v in histogram.items():
    #        thing.write("{}\t{}\n".format(k,v))
    #np.savetxt('timeseries.dat', np.column_stack((times, temps)))
    #### 

    return times, temps, w_temps, (nonw_temp_low, nonw_temp_high), (w_temp_low, w_temp_high)

    #temp_string = " ".join([str(temp) for temp in temps])
    #temp_string += " " + str(nonw_temp_low)

    #time_string = " ".join([str(time) for time in times])
    #time_string += " " + str(final_time)

    #if water_thermostat_style == 'proportional':
    #    w_temp_string = " ".join([str(temp) for temp in w_temps])
    #    w_temp_string += " " + str(w_temp_low)
    #    f.write("annealing-npoints\t={0} {0}\n".format(n_points+1))
    #    f.write("annealing-time\t={0} {0}\n".format(time_string))
    #    f.write("annealing-temp\t={0} {1}\n".format(temp_string, w_temp_string))
    #return possible_temps

def _generate_cooling_phase(time_start=25000, duration=50000, interval=5, 
        cooling_rate=500,
        nonw_temp_high=705, nonw_temp_low=305,
        w_temp_high=505, w_temp_low=305, water_thermostat_style=None):
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

    # This is a check to make sure water temperatures aren't 
    # excessively higher than the nonw temps
    if w_temp_high > nonw_temp_high:
        w_temp_high = nonw_temp_high
    if w_temp_low > nonw_temp_low:
        w_temp_low = nonw_temp_low

    if water_thermostat_style == 'proportional':
        scaling_factor = (w_temp_high - w_temp_low)/(nonw_temp_high - nonw_temp_low)
    w_temps = [w_temp_high]
    temps = [nonw_temp_high]
    times = [time_start]

    for i in np.arange(time_start, time_start+duration, interval):
        times.append(i)
        # If we've hit a cooling step, reduce the temperature range
        # Ensure that there is still a range of temperatures to choose from
        if i % cooling_rate == 0 and nonw_temp_high > nonw_temp_low+1:
            nonw_temp_high -= 1
        new_temp = np.random.randint(nonw_temp_low, nonw_temp_high)
        if water_thermostat_style == 'proportional':
            nonw_temp_change = new_temp - nonw_temp_low
            w_temp_change = scaling_factor * nonw_temp_change 
            new_w_temp = w_temp_low + w_temp_change
        elif water_thermostat_style =='plateau':
            if new_temp > w_temp_high:
                new_w_temp = w_temp_high
            else:
                new_w_temp  = new_temp
        else:
            sys.exit("Specify a proper water thermostatting style")

        temps.append(round(new_temp,2))
        w_temps.append(round(new_w_temp, 2))
        last_temp = new_temp
    return times, temps, w_temps, (nonw_temp_low, nonw_temp_high), (w_temp_low, w_temp_high)


def _write_annealing_lines(f, times, temps, w_temps):
    if len(times) == len(temps) == len(w_temps):
        temp_string = " ".join([str(temp) for temp in temps])
        w_temp_string = " ".join([str(temp) for temp in w_temps])
        time_string = " ".join([str(time) for time in times])
        f.write("annealing-npoints\t={0} {0}\n".format(len(times)))
        f.write("annealing-time\t={0} {0}\n".format(time_string))
        f.write("annealing-temp\t={0} {1}\n".format(temp_string, w_temp_string))
    else:
        sys.exit("Error, annealing lines do not have matching dimensions")





if __name__ == "__main__":
    main()

    
