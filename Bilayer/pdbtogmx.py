import math
import sys
import os
import numpy as np
import pdb

pdbfilename = 'dppccopy.pdb'
pdbfile = open(pdbfilename,'r')
pdblines = pdbfile.readlines()
outfile = open('dppccopy.gro','w')

for line in pdblines:
    stuff = line.split()
    print(stuff)
    atom_index = stuff[1]
    atom_name = stuff[2]
    res_name = stuff[3]
    resindex = float(stuff[4])
    xcoord = float(stuff[5])
    ycoord = float(stuff[6])
    zcoord = float(stuff[7])
    outfile.write("{:>5.0f}{:5s}{:>5s}{:>5s}{:8.3f}{:8.3f}{:8.3f}\n".format(resindex, res_name,
            atom_name, atom_index, xcoord/10, ycoord/10, zcoord/10))
            
