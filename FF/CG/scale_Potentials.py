import numpy as np
import os
import sys
import subprocess
import pdb
from optparse import OptionParser
""" Scale a tabulated potential"""

# These optimized pairs were derived in previous MSIBI work
optimized_pairs = [ ('P4', 'P4'), ('P4', 'P3'), ('P4', 'Nda'), ('P4', 'C1'),
                    ('P3', 'P3'), ('P3', 'Nda'), ('P3', 'C1'),
                    ('Nda', 'Nda'), ('Nda', 'C1'),
                    ('C1', 'C1') ] 


# Specify an input directory and an output directory
parser = OptionParser()
parser.add_option("--input", action = "store", type = "string", dest = "inp")
parser.add_option("--output", action = "store", type = "string", dest = "out")
parser.add_option("--scale", action = "store", type = "float", default = "1", dest = "scale")
parser.add_option("--ignore", action="store_true", default=False, dest = "ignore")
(options, args) = parser.parse_args()
if not options.inp or not options.out:
    sys.exit("ERROR: Specify input, option, and scale")


input_dir = str(os.getcwd())+ "/" + str(options.inp)
output_dir = str(os.getcwd())+  "/" + str(options.out)
scale = options.scale
current_dir = os.getcwd()

# Make directory
p = subprocess.Popen("mkdir -p {}".format(
            options.out), shell=True, stdout=subprocess.PIPE)
p.wait()

# Read every file in the input directory, 
dir_contents = os.listdir("{}".format(input_dir))
for pair_potential in dir_contents:
    pairs = pair_potential[:-4]
    first = pairs.split('-')[0].strip()
    second = pairs.split('-')[1].strip()
    input_file = os.path.join(input_dir , pair_potential)
    if options.ignore and (first,second) in optimized_pairs:
        p = subprocess.Popen("cp {} {}".format(input_file, output_dir), shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.wait()
        print(first,second)
        pass
    else: 
        tabulated_potential = np.loadtxt(input_file)
        # Scale energy and force terms
        tabulated_potential[:,1] *= scale
        tabulated_potential[:,2] *= scale
        # Print to output directory with same filename
        output_file = os.path.join(output_dir, pair_potential)
        np.savetxt(output_file, tabulated_potential)


