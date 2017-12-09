import math
import os
import random
random.seed(12345)

import numpy as np
from groupy.gbb import Gbb
from groupy.box import Box
from groupy.system import System
from groupy.builders.bilayer import Bilayer
import pdb

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", action="store", type="string", default="somebilayer", dest="filename")
parser.add_option("-a", "--APL", action="store",type="float", default=0.0, dest = "area_per_lipid")
parser.add_option("-r", "--rot", action="store", type ="float", default = 12.0, dest = "rotation")
parser.add_option("--DSPC", action="store",type="float", default=0.0, dest = "DSPC_frac")
parser.add_option("--DPPC", action="store",type="float", default=0.0, dest = "DPPC_frac")
parser.add_option("--acid12", action="store",type="float", default=0.0, dest = "acid12_frac")
parser.add_option("--acid16", action="store",type="float", default=0.0, dest = "acid16_frac")
parser.add_option("--acid22", action="store",type="float", default=0.0, dest = "acid22_frac")
parser.add_option("--alc12", action="store",type="float", default=0.0,  dest = "alc12_frac")
parser.add_option("--alc14", action="store",type="float", default=0.0, dest = "alc14_frac")
parser.add_option("--alc16", action="store",type="float", default=0.0, dest = "alc16_frac")
parser.add_option("--alc18", action="store",type="float", default=0.0, dest = "alc18_frac")
parser.add_option("--alc20", action="store",type="float", default=0.0, dest = "alc20_frac")
parser.add_option("--alc22", action="store",type="float", default=0.0, dest = "alc22_frac")
parser.add_option("--alc24", action="store",type="float", default=0.0, dest = "alc24_frac")
parser.add_option("--ISIS", action="store",type="float", default=0.0, dest = "isis_frac")
parser.add_option("--SS", action="store",type="float", default=0.0, dest = "ss_frac")
parser.add_option("--CHOL", action="store",type="float", default=0.0, dest = "chol_frac")
parser.add_option("--PMEA", action="store",type="float", default=0.0, dest = "pmea_frac")
parser.add_option("--water", action="store",type="float", default=0.0, dest = "water_frac")
(options, args) = parser.parse_args()

# load components
#base_path='/Users/yangah/Programs/setup/SimulationSetup_Remco/Prototypes/'
#base_path='Prototypes/'
base_path = '/raid6/homes/ahy3nz/Programs/setup/Bilayer/Prototypes'
dspc = Gbb(xml_prototype=base_path+'DSPC.xml',name='dspc')
dspc.rotate(angles=[0.0,0.0,-math.pi/4.0])
dppc = Gbb(xml_prototype=base_path+'DPPC.xml',name='dppc')
dppc.rotate(angles=[0.0,0.0,-math.pi/4.0])
acid12 = Gbb(xml_prototype=base_path+'c12-acid.xml',name='acid12')
acid16 = Gbb(xml_prototype=base_path+'c16-acid.xml',name='acid16')
acid22 = Gbb(xml_prototype=base_path+'c22-acid.xml',name='acid22')
alc12 = Gbb(xml_prototype=base_path+'c12-alcohol.xml',name='alc12')
alc14 = Gbb(xml_prototype=base_path+'c14-alcohol.xml',name='alc14')
alc16 = Gbb(xml_prototype=base_path+'c16-alcohol.xml',name='alc16')
alc18 = Gbb(xml_prototype=base_path+'c18-alcohol.xml',name='alc18')
alc20 = Gbb(xml_prototype=base_path+'c20-alcohol.xml',name='alc20')
alc22 = Gbb(xml_prototype=base_path+'c22-alcohol.xml',name='alc22')
alc24 = Gbb(xml_prototype=base_path+'c24-alcohol.xml',name='alc24')
isis = Gbb(xml_prototype=base_path+'isostearylisostearate.xml',name='isis') 
isis.rotate(angles=[math.pi,0.0,0.0])
ss = Gbb(xml_prototype=base_path+'ss.xml',name='ss') 
ss.rotate(angles=[math.pi,0.0,0.0])
chol = Gbb(xml_prototype=base_path+'CHOL.xml',name='chol')
chol.rotate(angles=[math.pi,0.0,0.0])
pmea = Gbb(xml_prototype=base_path+'PMEA.xml',name='pmea')
pmea.rotate(angles=[0.0,-math.pi/2,0.0])
water = Gbb(xml_prototype=base_path+'spc.xml',name='water')
#DSPC thickness really 16
lipids = [(dspc,  options.DSPC_frac , 16.0),
          (dppc,  options.DPPC_frac , 13.0),
          (acid12,options.acid12_frac , 13.0),
          (acid16,options.acid16_frac , 22.0),
          (acid22,options.acid22_frac , 16.0),
          (alc12, options.alc12_frac , 13.0),
          (alc14, options.alc14_frac , 14.0),
          (alc16, options.alc16_frac , 13.0),
          (alc18, options.alc18_frac , 12.0),
          (alc20, options.alc20_frac , 15.0),
          (alc22, options.alc22_frac , 16.0),
          (alc24, options.alc24_frac , 16.0),
          (isis,  options.isis_frac , 15.0),
          (ss,    options.ss_frac , 15.0),
          (chol,  options.chol_frac , 20.0),
          (pmea,  options.pmea_frac , 22.0),
          (water, options.water_frac , 0.0)]
'''
lipids = [(dspc,  options.DSPC_frac , 25.0),
          (dppc,  options.DPPC_frac , 13.0),
          (acid16,options.acid16_frac , 25.0),
          (acid22,options.acid22_frac , 14.0),
          (alc12, options.alc12_frac , 16.0),
          (alc14, options.alc14_frac , 17.0),
          (alc16, options.alc16_frac , 16.0),
          (alc18, options.alc18_frac , 15.0),
          (alc20, options.alc20_frac , 18.0),
          (alc22, options.alc22_frac , 19.0),
          (alc24, options.alc24_frac , 19.0),
          (isis,  options.isis_frac , 18.0),
          (ss,    options.ss_frac , 18.0),
          (chol,  options.chol_frac , 23.0),
          (pmea,  options.pmea_frac , 25.0),
          (water, options.water_frac , 0.0)]
          '''

lipid_positions = [range(18), [], [], [], range(18, 36), [], [], [], [], [], []]    # CHANGE THIS!!!

phasesep = 'no'
#max_rot = 24.0 * math.pi / 180 * 2
#Try no rotation here 
#max_rot =6 * math.pi / 180 * 2
max_rot = options.rotation * math.pi/180

n_x = n_y = 8
solvent_per_lipid=20
random_z_displacement = 3.0
#area per lipid was initially 44 for remco
if (phasesep=='no'):
    print('No phase separation')
    bilayer = Bilayer(lipids, n_x=n_x, n_y=n_y,
        area_per_lipid=options.area_per_lipid, mirror=False,
        solvent_per_lipid=20, solvent=water)

    
    resid_atom_list = [ (dspc, 2 * n_x * n_y * np.float(lipids[0][1])),
                        (dppc, 2 * n_x * n_y * np.float(lipids[1][1])),
                        (acid12, 2 * n_x * n_y * np.float(lipids[2][1])),
                        (acid16, 2 * n_x * n_y * np.float(lipids[3][1])),
                        (acid22, 2 * n_x * n_y * np.float(lipids[4][1])),
                        (alc12, 2 * n_x * n_y * np.float(lipids[5][1])),
                        (alc14, 2 * n_x * n_y * np.float(lipids[6][1])),
                        (alc16, 2 * n_x * n_y * np.float(lipids[7][1])),
                        (alc18, 2 * n_x * n_y * np.float(lipids[8][1])),
                        (alc20, 2 * n_x * n_y * np.float(lipids[9][1])),
                        (alc22, 2 * n_x * n_y * np.float(lipids[10][1])),
                        (alc24, 2 * n_x * n_y * np.float(lipids[11][1])),
                        (isis, 2 * n_x * n_y * np.float(lipids[12][1])),
                        (ss, 2 * n_x * n_y * np.float(lipids[13][1])),
                        (chol, 2 * n_x * n_y * np.float(lipids[14][1])),
                        (pmea, 2 * n_x * n_y * np.float(lipids[15][1])),
                        (water, 2 * n_x * n_y * np.float(solvent_per_lipid)) ]


    for lipid in bilayer.molecules[:n_x*n_y*2]:
        rot_z = random.random() * math.pi * 2
        lipid.rotate(angles=[0.0, 0.0, rot_z], fixed_atom=12)
        rot_x = 0.0
        rot_y = max_rot
        rot_z = 0.0 #random.random() * math.pi * 2
        lipid.rotate(angles=[rot_x, rot_y, rot_z], fixed_atom=12)
        lipid.calc_com()
        disp_z = random.random() * random_z_displacement * -np.sign(lipid.com[2])
        lipid.translate([0, 0, disp_z])
#    system = System(box=bilayer.box, gbbs=bilayer.molecules)
    system=System(box=bilayer.box, gbbs=bilayer.molecules)
    system.set_resids(bilayer.molecules)
    #need to make a list for the names of bilayer molecules
    resid_name_list=list()
    for i in bilayer.molecules:
        resid_name_list.append(i.name)
    system.set_resnames(resid_name_list)

    #mygbb = Gbb() 
    #mygbb.load_lammps_data(data_file='bilayer.lammpsdata')
    #pdb.set_trace()
    #system.sort_by_name([x[0].name for x in lipids])
    #Break, check what system has, print to gromacs not lammps, need
    #to add functionality in system.py
    system.print_lammpsdata(filename=(options.filename+'.lammpsdata'),
       ff_param_set='gromos53a6')
    #mygbb = Gbb() 
    #mygbb.load_lammps_data(data_file='bilayer.lammpsdata')
    #system.print_grodata(resid_atom_list, mygbb=mygbb, filename='bilayer.gro', ff_param_set='gromos53a6')
    #os.system('vmd_MACOSXX86 -e vis.vmd')
else:
    print('Phase separation')
    bilayer = Bilayer(lipids, n_x=6, n_y=6,
        area_per_lipid=45, mirror=True,
        solvent_per_lipid=10, solvent=water, lipid_positions=lipid_positions)
    for lipid in bilayer.molecules[36:72]:
        lipid.rotate([0.0, 0.0, np.pi])
    system = System(box=bilayer.box, gbbs=bilayer.molecules)
    system.sort_by_name([x[0].name for x in lipids])
    #system.sort_by_name([dppc.name,dspc.name,acid16.name,acid22.name,alc16.name,alc22.name,alc18.name,isis.name,chol.name,water.name])
    system.print_lammpsdata(filename='bilayer.lammpsdata',
            ff_param_set='gromos53a6')
#    os.system('vmd_MACOSXX86 -e vis.vmd')

