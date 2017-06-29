import numpy as np
import os
import sys
import subprocess
import pdb
from optparse import OptionParser
""" Scale a tabulated potential"""


# Specify an input directory and an output directory
parser = OptionParser()
parser.add_option("--input", action = "store", type = "string", dest = "inp")
parser.add_option("--output", action = "store", type = "string", dest = "out")
parser.add_option("--scale", action = "store", type = "float", default = "1", dest = "scale")
(options, args) = parser.parse_args()
#if not options.inp or not options.out or not options.scale:
#    sys.exit("ERROR: Specify input, option, and scale")


input_dir = str(os.getcwd())+ "/" + str(options.inp)
output_dir = str(os.getcwd())+  "/" + str(options.out)
scale = options.scale
current_dir = os.getcwd()

# Make directory
p=subprocess.Popen("mkdir -p {}".format(
            options.out), shell=True, stdout=subprocess.PIPE)
p.wait()

# Read every file in the input directory, 
dir_contents = os.listdir("{}".format(input_dir))
for pair_potential in dir_contents:
    input_file = os.path.join(input_dir , pair_potential)
    tabulated_potential = np.loadtxt(input_file)
    # Scale energy and force terms
    tabulated_potential[:,1] *= scale
    tabulated_potential[:,2] *= scale
    # Print to output directory with same filename
    output_file = os.path.join(output_dir, pair_potential)
    np.savetxt(output_file, tabulated_potential)


