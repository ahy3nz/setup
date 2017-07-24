import os
import pdb
import numpy as np
import sys
import subprocess

filename = 'martini_pair.txt'
#filename='testpair.txt'
martini_file = open(filename,'r')
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
        subprocess.Popen("python generate_Hoomd_Table.py -f {}-{} --c6 {} --c12 {}".format(
            atom_1, atom_2, c6, c12), shell=True, stdout=subprocess.PIPE)
