import math
import sys 
import os
import numpy as np
import pdb
from optparse import OptionParser

'''
The goal of this is to fill a gro file with the coordinates of a respective lmps file
A priori coordinates in the gro file are irrelevant, 
but still need correct atom and indices. Topology file remains untouched
'''
parser = OptionParser()
parser.add_option("--gro", action = "store", type = "string", dest = "grofilename")
parser.add_option("--lmps", action = "store", type = "string", dest = "lmpsfilename")
(options, args) = parser.parse_args()

grofilename = options.grofilename
lmpsfilename = options.lmpsfilename
#grofilename = 'LammpsDSPC-100.gro'
#lmpsfilename = 'LammpsDSPC-100.lammpsdata'

grofile = open(str(grofilename),'r')
lmpsfile = open(str(lmpsfilename),'r')

grolines = grofile.readlines()
lmpslines = lmpsfile.readlines()

n_atoms = int(grolines[1])
print('There are {} atoms'.format(n_atoms))
lmps_xyz = [None] * (n_atoms+1)
gro_xyz = [list()] * (n_atoms+1)
#gro_xyz=list()
lmps_atom_list_start = 0 
#Get Box stuff
if "zlo" in lmpslines[16]:
    lmps_z_box = float(lmpslines[16].split()[1])
    lmps_y_box = float(lmpslines[15].split()[1])
    lmps_x_box = float(lmpslines[14].split()[1])
        
    #print (lmpslines[15])
else:
    print('Problem reading box in lammps')
    sys.exit()

gro_z_box = lmps_z_box/5
gro_y_box = lmps_y_box/5
gro_x_box = lmps_x_box/5
#Read lmps lines, looking for the section that contains the atom coordiantes
for i, line in enumerate(lmpslines):
    #if "Atoms" in line.split()[0]:
    if "Atoms" in line:
        lmps_atom_list_start = i+1
    else:
        pass

#Start looking at coordinates, scaling and transforming for gro coordinates
#Index of gro_xyz corresponds to Lammps indices, starting at index 1
for i in range(lmps_atom_list_start+1, lmps_atom_list_start + n_atoms + 1):
    if len(lmpslines[i]) == 0:
        pass
    else:
        (x,y,z) = ( float(lmpslines[i].split()[4])/10,
                    float(lmpslines[i].split()[5])/10,
                    (float(lmpslines[i].split()[6]) + lmps_z_box)/10)
        #Lmps coordinates start at index 1
        atom_index = int(lmpslines[i].split()[0])
        gro_xyz[atom_index] = (x,y,z)

for i, line in enumerate(grolines):
    if i == 0 or i == 1:
        pass
    elif i < n_atoms+2:
        splitline = line.split()
        firstpart = line[:20]
        atom_index = i - 1
        splitline[-1] = gro_xyz[atom_index][2]
        splitline[-2] = gro_xyz[atom_index][1]
        splitline[-3] = gro_xyz[atom_index][0]
        #grolines[i] = splitline
        grolines[i] = '{}{:8.3f}{:8.3f}{:8.3f}'.format(firstpart,
                gro_xyz[atom_index][0],
                gro_xyz[atom_index][1],
                gro_xyz[atom_index][2],)


grolines[-1] = '{:10.6f} {:10.6f} {:10.6f}'.format(gro_x_box, gro_y_box, gro_z_box)
print("Writing to DidItWork.gro")

outfile = open('DidItWork.gro','w')
for i,line in enumerate(grolines):
    outfile.write(grolines[i]+'\n')
    #print(line)
    #if i <= 1:
    #    outfile.write(''.join(str(x) for x in line))
    #elif len(line) ==3:
    #    outfile.write('{:10.5f}{:10.5f}{:10.5f}\n'.format(line[0],line[1],line[2]))
    #elif i > 1:
    #    outfile.write('{:>5.0f}{:5s}{:>5s}{:>5s}{:8.3f}{:8.3f}{:8.3f}\n'.format(
    #        float(line[0]), line[1], line[2], line[3], line[4],line[5] ))
    #else:
    #    outfile.write(''.join(str(x)for x in line)+'\n')
