#import xml.etree.ElementTree as ET
from collections import OrderedDict
from lxml import etree as ET
import itertools as it

#beadtypes = ['W', 'HS', 'E1','C2', 'C3', 'PCP', 'PCN', 'OH5' , 'head']
beadtypes = OrderedDict()
beadtypes['W'] = '72'
beadtypes['HS'] = '36'
beadtypes['E1'] = '58.04'
beadtypes['C2'] = '29.06'
beadtypes['C3'] = '42.08'
beadtypes['PCP'] =  '123.03'
beadtypes['PCN'] =  '73.14'
beadtypes['OH5'] =  '17.01'
beadtypes['head'] = '44.00'
root = ET.Element('ForceField')
atomtypes = ET.SubElement(root, "AtomTypes")
for bead, mass in beadtypes.items():
    new_node = ET.SubElement(atomtypes, "Type")
    new_node.attrib['name'] = bead
    new_node.attrib['class'] = bead
    new_node.attrib['element'] = "_"+bead
#    if bead == 'HS':
#        new_node.attrib['mass'] = "36"
#    else:
#        new_node.attrib['mass'] = "72"
    new_node.attrib['mass'] = mass
    new_node.attrib['def'] = '[_{}]'.format(bead)
    atomtypes.append(new_node)

bonds = ET.SubElement(root, "HarmonicBondForce")
for a, b in it.combinations_with_replacement(beadtypes, 2):
    new_node = ET.SubElement(bonds, "Bond")
    new_node.attrib['class1'] = a
    new_node.attrib['class2'] = b
    new_node.attrib['length'] = "1"
    new_node.attrib['k'] = "1"

angles = ET.SubElement(root, "HarmonicAngleForce")
for b in beadtypes:
    for a, c in it.combinations_with_replacement(beadtypes, 2):
        new_node = ET.SubElement(angles, "Angle")
        new_node.attrib['class1'] = a
        new_node.attrib['class2'] = b
        new_node.attrib['class3'] = c
        new_node.attrib['angle'] = "1"
        new_node.attrib['k'] = "1"

nb = ET.SubElement(root, "NonbondedForce")
nb.attrib['coulomb14scale'] = '0'
nb.attrib['lj14scale'] = '0'
for bead in beadtypes:
    new_node = ET.SubElement(nb, "Atom")
    new_node.attrib['type'] = bead
    new_node.attrib['charge'] = "0"
    new_node.attrib['sigma'] = "1"
    new_node.attrib['epsilon'] = "1"

tree = ET.ElementTree(root)

tree.write("out.xml", pretty_print=True)
