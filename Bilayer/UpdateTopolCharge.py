import math
import sys 
import os
import numpy as np
import pdb
from optparse import OptionParser

'''
Change the #include lines of the topology
to go from ff_b.itp to ff.itp
'''
parser = OptionParser()
parser.add_option("-p", action = "store", type = "string", dest = "topfilename")
(options, args) = parser.parse_args()
topfile = open(options.topfilename,'r')
topfilelines = topfile.readlines()
newlines = []
for i, line in enumerate(topfilelines):
    if "ff_b.itp" in line or "spc_b.itp" in line:
        newline = (";" + line)
    elif "ff.itp" in line or "spc.itp" in line:
        newline = line.replace(";", "")
    else:
        newline = line
    newlines.append(newline)
topfile.close()
topfile = open(options.topfilename,'w')
for i, line in enumerate(newlines):
    topfile.write(line)


