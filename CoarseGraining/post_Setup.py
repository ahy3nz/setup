import pdb
import numpy as np
import mdtraj as mdtraj
import sys
from Prototypes_CG import *
from optparse import OptionParser


"""
Due to lack of atomtyping in foyer for these CG systems,
go back and rewrite .gro files, accounting for 
residue naming and indexing. 
Requires table of contents file that lists
(residue index, residue name, n_atoms in that residue)
"""

def split_lines(source = None):
    """ Simple script to take a file.readlines() and split each line into a list

    Parameters
    ----------
    source : list
        A file whose lines have been read into a list

    Returns
    -------
    destination : list
        A list whose entries are lists of the split lines
    """
    destination = []
    for i, line in enumerate(source):
        destination.append(line.split())
    return destination

def gather_header(old_gro_lines = None, new_gro_file = None):
    """ Write gro file header information

    Parameters
    ---------
    old_gro_lines : list
        List whose entries are split lines of gro file
    new_gro_file : file
        Output gro file

    Returns
    -------
    old_gro_lines : list
        List whose entries are split lines of gro file, with updated lines removed

        """
    line_index = 0
    while (len(old_gro_lines[0]) != 6 or line_index < 2):
        new_gro_file.write(" ".join(old_gro_lines.pop(0)) + "\n")
        line_index += 1
    return old_gro_lines

def gather_body(old_gro_lines = None, table_of_contents = None, new_gro_file = None):
    """ Write body of the topology file

    Parameters
    ----------
    old_gro_lines : list
        List whose entries are split lines of gro file
    table_of_contents : file
        Input file listing residue index, residue name, n_particles for tha residue
    new_gro_file : file
        Output gro file

    Returns
    -------
    old_gro_lines : list
        List whose entries are split lines of gro file, with updated lines removed

    Notes
    -----
    Gro file format
    {:>5.0f}{:5s}{:>5s}{:>5s}{:8.3f}{:8.3f}{:8.3f}
    resindex(1 indexed), resname, atomtype, atomindex(1 indexed), x-coord, y-coord, z-coord

    """
    # Iterate through the table of contents, which then tells the inner loop how many times to iterate
    for i in range(len(table_of_contents)-1):
        entry = table_of_contents[i]
        res_index = entry[0]
        res_name = entry[1]
        res_atoms = int(entry[2])
        for j in range(res_atoms):
            old_line = old_gro_lines.pop(0)
            new_line = (float(res_index), res_name, old_line[1], old_line[2], float(old_line[3]),
                    float(old_line[4]), float(old_line[5]))
            new_gro_file.write("{:>5.0f}{:5s}{:>5s}{:>5s}{:8.3f}{:8.3f}{:8.3f}\n".format(new_line[0], new_line[1],
                new_line[2], new_line[3], new_line[4], new_line[5], new_line[6]))
    return old_gro_lines

def gather_footer(table_of_contents= None, new_gro_file = None):
    """ Print footer information (box size)

    Parameters
    ----------
    table_of_contents : file
        Input file listing residue index, residue name, n_particles for tha residue
    new_gro_file : file
        Output gro file

    """
    new_gro_file.write(' '.join(table_of_contents[-1])+"\n")

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default = "CG_bilayer", dest = "filename")
(options, args) = parser.parse_args()

filename = options.filename
table_of_contents_file = open(filename+'.dat', 'r').readlines()
old_gro_file = open(filename + '.gro', 'r').readlines()
new_gro_file = open(filename + '.gro', 'w')
print("Updating <{0}.gro> ...".format(filename))

# Split each line in the old gro file
old_gro_lines = split_lines(old_gro_file)

# Split each line in the table of contents file
table_of_contents = split_lines(table_of_contents_file)

# Header information
old_gro_lines = gather_header(old_gro_lines = old_gro_lines, new_gro_file = new_gro_file)
        
# Body information 
old_gro_lines = gather_body(old_gro_lines = old_gro_lines, 
        table_of_contents = table_of_contents, new_gro_file = new_gro_file)

# Footer information
gather_footer(table_of_contents = table_of_contents, new_gro_file = new_gro_file)

