import os
import pdb
import numpy as np
import sys
import subprocess
import argparse
import ff_utils

""" Generate tabulated potnetials from a pair file
First take c6 and c12 to get sigma and epsilon
Scale sigma and epsilon if specified
Generate table (morse or LJ)
"""
# Optimized pairs are those derived in earlier work by Tim
optimized_pairs = [ ('P4', 'P4'), ('P4', 'P3'), ('P4', 'Nda'), ('P4', 'C1'),
                    ('P3', 'P3'), ('P3', 'Nda'), ('P3', 'C1'),
                    ('Nda', 'Nda'), ('Nda', 'C1'),
                    ('C1', 'C1') ] 


parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename", help="Filename", default="martini_pair.dat")
parser.add_argument("--scale_sigma", type=float,help="Scale sigma val", default=1.0)
parser.add_argument("--scale_epsilon", type=float,help="Scale epsilon val", 
        default=1.0)
parser.add_argument("--output", default="new_pairs")
parser.add_argument("--plot", action="store_true", default=False)
parser.add_argument("--morse", action="store_true", default=False)
parser.add_argument("--lj", action="store_true", default=False)
args = parser.parse_args()

p = subprocess.Popen("mkdir -p {}".format(args.output),shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p.wait()

#filename = 'martini_pair.dat'
#filename='testpair.txt'
if not args.lj and not args.morse:
    sys.exit("Specify LJ and/or morse potnetial fits")
martini_file = open(args.filename,'r')
martini_lines = martini_file.readlines()
for i, line in enumerate(martini_lines):
    split_lines = line.split(';')
    split_lines = split_lines[0].split()
    # element 0 is first atom
    if(line.rstrip()):
        print(split_lines)
        atom_1 = split_lines[0]
        # element 1 is second atom
        atom_2 = split_lines[1]
        # elemtn 3 is  c6 term
        c6 = split_lines[3]
        # element 4 is c12term
        c12 = split_lines[4]

        # For computing sig eps:
        sig, eps = ff_utils.gmx_to_sigeps(c6=float(c6), c12=float(c12), type_A=atom_1, type_B=atom_2)
        scaled_sig = sig*args.scale_sigma
        scaled_eps = eps*args.scale_epsilon
        new_c6 = 4*scaled_eps*(scaled_sig**6)
        new_c12 = 4*scaled_eps*(scaled_sig**12)

        
        print("Sig {} -> {}".format(round(sig,3), round(scaled_sig,3)))
        print("Eps {} -> {}".format(round(eps,3), round(scaled_eps,3)))
        # For making tables:
        if args.lj:
            ff_utils.generate_Table(c6=new_c6, c12=new_c12, 
                output='{}/{}-{}'.format(args.output, atom_1, atom_2),plot=args.plot)
        if args.morse:
            # Get some morse potentials
            morse_params = ff_utils.LJ_to_morse(start_fit=0.9*scaled_sig, end_fit=1.1*scaled_sig,
                sigma=scaled_sig, eps=scaled_eps, output_plot=args.plot, output_name="{}/{}-{}_fit_normal".format(args.output,atom_1, atom_2))

            ff_utils.generate_morse_Table(**morse_params,
                output='{}/{}-{}'.format(args.output, atom_1, atom_2))

        
