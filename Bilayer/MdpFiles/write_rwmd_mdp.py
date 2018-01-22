import os
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
        time_anneal=25000, final_time=50000):
    """ Generate annealing points
    temp_low : lowest temperature, float
    temp_high : highest temperature, float
    interval : duration for each annealing interval, int (ps)
    time_anneal: duration of annealing, int (ps) 
    """
    n_points = int(time_anneal/interval)
    times = [0]
    temps = [temp_low]
    last_temp = 305
    for i in np.arange(1,n_points):
        # If we're at the lowest temperature, only go up
        if abs(last_temp - temp_low) < 0.1:
            dtemp = 5 * np.random.randint(0,2)
        # If we're at the highest temperature, only go down
        elif abs(last_temp - temp_high) < 0.1:
            dtemp = 5 * np.random.randint(-1,1)
        else:
            dtemp = 5 * np.random.randint(-1,2)
        time = i*interval
        new_temp = last_temp + dtemp
        times.append(time)
        temps.append(new_temp)
        last_temp = new_temp

    temp_string = " ".join([str(temp) for temp in temps])
    temp_string += " " + str(temp_low)
    time_string = " ".join([str(time) for time in times])
    time_string += " " + str(final_time)
    f.write("annealing-npoints\t={0} 2 \n".format(n_points+1))
    f.write("annealing-time\t={0} 0 {1}\n".format(time_string, str(final_time)))
    f.write("annealing-temp\t={0} {1} {1}\n".format(temp_string, str(temp_low)))

if __name__ == "__main__":
    with open('RW.mdp','w') as f:
        _write_body(f, ref_temp=305)
        _generate_annealing_points(f, temp_low=305, temp_high=560, interval=5,
                time_anneal=25000, final_time=50000)
