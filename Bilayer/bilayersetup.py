import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "somebilayer", dest = "filename")
parser.add_option("-a", "--APL", action="store",type="float", default = 0.0, dest = "area_per_lipid")
parser.add_option("-r", "--rot", action="store", type ="float", default = 12.0, dest = "rotation")
parser.add_option("--max", action="store", type="float",  dest = maxtemp)
parser.add_option("--min", action="store", type="float",  dest = mintemp)
parser.add_option("--DSPC", action="store",type="float", default = 0.0, dest = "DSPC_frac")
parser.add_option("--DPPC", action="store",type="float", default = 0.0, dest = "DPPC_frac")
parser.add_option("--acd16", action="store",type="float", default = 0.0, dest = "acid16_frac")
parser.add_option("--acd22", action="store",type="float", default = 0.0, dest = "acid22_frac")
parser.add_option("--alc12", action="store",type="float", default = 0.0,  dest = "alc12_frac")
parser.add_option("--alc14", action="store",type="float", default = 0.0, dest = "alc14_frac")
parser.add_option("--alc16", action="store",type="float", default = 0.0, dest = "alc16_frac")
parser.add_option("--alc18", action="store",type="float", default = 0.0, dest = "alc18_frac")
parser.add_option("--alc20", action="store",type="float", default = 0.0, dest = "alc20_frac")
parser.add_option("--alc22", action="store",type="float", default = 0.0, dest = "alc22_frac")
parser.add_option("--alc24", action="store",type="float", default = 0.0, dest = "alc24_frac")
parser.add_option("--ISIS", action="store",type="float", default = 0.0, dest = "isis_frac")
parser.add_option("--SS", action="store",type="float", default = 0.0, dest = "ss_frac")
parser.add_option("--CHOL", action="store",type="float", default = 0.0, dest = "chol_frac")
parser.add_option("--PMEA", action="store",type="float", default = 0.0, dest = "pmea_frac")
parser.add_option("--water", action="store",type="float", default = 0.0, dest = "water_frac")
(options, args) = parser.parse_args()

# Write out initial parameters
outfile = open((options.filename + 'initparam.txt'),'w')
outfile.write('Initial APL: {}\n'.format(options.area_per_lipid))
outfile.write('Initial Tilt: {}\n'.format(options.rotation))
outfile.close()

# Write Cori scripts
os.system("python WriteCoriScript.py -f {}".format(options.filename))

# Write the system in lammps
os.system(("python init-bilayer_tilted.py -f {} -a {} -r {} --DSPC {} --DPPC {} --acid16 {} "
            "--acid22 {} --alc12 {} --alc14 {} --alc16 {} --alc18 {} --alc20 {} --alc22 {} "
            "--alc24 {} --ISIS {} --SS {} --CHOL {} --PMEA {} --water {}").format( options.filename,
                                options.area_per_lipid,
                                options.rotation,
                                options.DSPC_frac,
                                options.DPPC_frac,
                                options.acid16_frac,
                                options.acid22_frac,
                                options.alc12_frac,
                                options.alc14_frac,
                                options.alc16_frac,
                                options.alc18_frac,
                                options.alc20_frac,
                                options.alc22_frac,
                                options.alc24_frac,
                                options.isis_frac,
                                options.ss_frac,
                                options.chol_frac,
                                options.pmea_frac,
                                options.water_frac))


# Write the system in gromacs 
os.system(("python bilayer_lmps2gmx.py -f {} -a {} -r {} --DSPC {} --DPPC {} --acid16 {} "
            "--acid22 {} --alc12 {} --alc14 {} --alc16 {} --alc18 {} --alc20 {} --alc22 {} "
            "--alc24 {} --ISIS {} --SS {} --CHOL {} --PMEA {} --water {}").format( options.filename,
                             options.area_per_lipid,
                             options.rotation,
                             options.DSPC_frac,
                             options.DPPC_frac,
                             options.acid16_frac,
                             options.acid22_frac,
                             options.alc12_frac,
                             options.alc14_frac,
                             options.alc16_frac,
                             options.alc18_frac,
                             options.alc20_frac,
                             options.alc22_frac,
                             options.alc24_frac,
                             options.isis_frac,
                             options.ss_frac,
                             options.chol_frac,
                             options.pmea_frac,
                             options.water_frac))

# Write Simulated Tempering mdp file
os.system("python writeSTmdp.py -f {} --min {} --max {}".format(options.filename, options.mintemp, options.maxtemp))

# Run equilibration steps
os.system("bash mdpreppbs.sh {}".format(options.filename))

# Sort files into folders
os.system('mkdir -p {}'.format(options.filename))
os.system('mv *{}* {}'.format(options.filename, options.filename))


