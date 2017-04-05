import itertools
import os
import sys
import numpy as np

atomtypes = ['P5' , 'P4' , 'BP4' , 'SP4' , 'P3'	, 'P2', 'P1','SP1' ,
             'Nda' 	,'Nd', 'SNd', 'Na','SNa', 'N0'	,
             'C5', 'SC5', 'C4', 'SC4', 'C3'	, 'SC3'	, 'C2', 'AC2', 'SC2', 'C1', 'AC1', 'SC1',
             'Qda', 'Qd', 'Qa', 'Q0']
outfile = open('xmlstuff.xml','w')
outfile.write("<Forcefield>\n")

# Atomtype xml writing
outfile.write("<Atomtypes>\n")
for x in atomtypes:
    # Exception is Q0(charge +1.0) and Qa(charge -1.0)
    if "Q0" == x:
        xml_line="\t<Type name=\"Q0\" class=\"Q0\" element=\"_Q0\" mass=\"72\" charge=\"1.0\" def=\"[_Q0]\" desc=\"PC head\"/>"
    elif "Qa" == x:
        xml_line="\t<Type name=\"Qa\" class=\"Qa\" element=\"_Qa\" mass=\"72\" charge=\"-1.0\" def=\"[_Qa]\" desc=\"PC head\"/>"
    else:
        xml_line = "\t<Type name=\"{0}\" class=\"{0}\" element=\"_{0}\" mass=\"72\" def=\"[_{0}]\"/>".format(x)
    outfile.write(xml_line+"\n")
outfile.write("</Atomtypes>\n")

# Bond xml writing
outfile.write("<HarmonicBondForce>\n")
for x,y in itertools.combinations_with_replacement(atomtypes,2):
    # exception bond, length 2 and k2, this is for Na-Na in a PC molecule
    if 'Na' == x and 'Na' == y:
        xml_line ="\t<Bond class1=\"{}\" class2=\"{}\" length=\"2\" k=\"2\"/>".format(x,y)
    # Standard bond, length 1 and k 1
    else:
        xml_line = "\t<Bond class1=\"{}\" class2=\"{}\" length=\"1\" k=\"1\"/>".format(x,y)
    outfile.write(xml_line+"\n")
outfile.write("</HarmonicBondForce>\n")

#<Bond class1="tail" class2="tail" length="1" k="1"/>
outfile.write("<HarmonicAngleForce>\n")
# Angle xml writing
# 1 -2 -3 is the same as 3-2-1, 1-3-2 is different
# The middle atom is very crucial, but the outside two order doens't matter
# Could do combination with replacement of the outer two atoms, manually iterate through all potnetial inner atoms
# Q0-Na-Na the exception angle
for x,z in itertools.combinations_with_replacement(atomtypes,2):
    for y in atomtypes: 
        if ('Na'== x and 'Na' == y and 'Q0' == z) or ('Q0' == x and 'Na' == y and 'Na' == z):
            xml_line = "\t<Angle class1=\"{}\" class2=\"{}\" class3=\"{}\" angle=\"2\" k=\"2\"/>".format(x,y,z)
        else:
            xml_line = "\t<Angle class1=\"{}\" class2=\"{}\" class3=\"{}\" angle=\"1\" k=\"1\"/>".format(x,y,z)
        outfile.write(xml_line+"\n")
#<Angle class1="tail" class2="mhead2" class3="oh2" angle="7" k="7"/>
outfile.write("</HarmonicAngleForce>\n")
outfile.write("</Forcefield\n")
outfile.close()

