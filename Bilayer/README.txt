init-bilayer_tilted.py
Constructs bilayer systems in lammps while using option parser parameters, rotations still need specification

bilayer_lmps2gmx.py
Constructs bilayer systmes in gmx while using option parser parameters. Much of the same construction as the lammps construction

bilayersetup.py
Calls init-bilayer_tilted and bilayer_lmps2gmx.py using option parser parameters, writes out initial structure
parameters, writes slurm submission script to cori, casll mdprep.sh (see below) to prepare for simulated
tempering or classical MD

Option parser parameters
-f : specify a filename for lammps file, gmx gro file, and gmx top file 
-a : specify an area per lipid for bilayer construction
--DSPC : specify the DSPC fraction (1 for 100%, 0.33 for 33%, 0.5 for 50%, etc)
--alc12: specify the alcohol fraction
--acid16: specify the acid fraction
...etc for all other major components

mdp files
These mdp files are used for the various steps in gromacs simulations
energy minimization, NVT equilibration, and NPT production runs
The general workflow is outlined in mdprep.sh

mdprep.sh
Call this shell script, entering the filename as the first argument
Runs grompp for EM, mdruns for EM, grompps for NVT equilibration, mdruns for NVT eq, and grompps for NPT production
Ex: bash mdprep.sh bilayer (this will run md prep for the gro and top file called "bilayer")

Note regarding top files:
The top files have a "#include" statement that points to a file called ff.itp
These top files need to be edited so that the directories in the include statement point to the absolute directory of this ff.itp file
The ff.itp file is located in ./gromos53a6/ff.itp, but absolute directory needs to be listed in the topology file
The "#include gromos53a6.ff/forcefield.itp" should be commented out
The ff.itp points to the topologies of some common prototypes for gel-phase bilayers

common_apl.txt
This is a list of some possible APL to initialize the bilayers at. If APLs are not correctly specified, the EM step will likely crash due to
atoms on top of each other

UpdateGroWithLmps.py
This is a supplemental file to replace the coordinates of a gro file with the coordinates of a lammps file. The topology file is untouched.
The only requirement is that the atoms and types of the gro file match the atoms and types of the lammps file so the gro file can easily
have respective coordinates replaced. Arguments are --gro and --lmps for the respective gro and lammps files

UpdateTopolCharge.py
the mdprep framework uses a slowgrowth method for energy minimization on the vdw forces, this file changes the
topology files to include charges for equilibration post-EM
