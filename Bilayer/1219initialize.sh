#!/bin/bash
# (1:1) DSPC-C12OH
python bilayersetup.py -f DSPC-50_alc12-50_1219 --DSPC 0.50 --alc22 0.50 -a 31 -r 9

# (1:4) DSPC-C24OH
python bilayersetup.py -f DSPC-20_alc24-80_1219 --DSPC 0.20 --alc24 0.80 -a 26 -r 20

# (1:2) DSPC-C12OH
python bilayersetup.py -f DSPC-34_alc12-66_1219 --DSPC 0.34 --alc12 0.66 -a 29 -r 18

# (1:2) DSPC-C22FFA
python bilayersetup.py -f DSPC-34_acd22-66_1219 --DSPC 0.34 --acid22 0.66 -a 29 -r 18

# (1:1:1) DSPC-C22OH-C22FFA
python bilayersetup.py -f DSPC-34_alc22-33_acd22-33_1219 --DSPC 0.34 --alc22 0.33 --acid22 0.33 -a 24 -r 6

# (1:1:1) DSPC-C24OH-C16FFA
python bilayersetup.py -f DSPC-34_alc24-33_acd16-33_1219 --DSPC 0.34 --alc24 0.33 --acid16 0.33 -a 28 -r 14

# (1:2) DSPC-C24OH
python bilayersetup.py -f DSPC-34_alc24-66_1219 --DSPC 0.34 --alc24 0.66 -a 28 -r 14

# (1:2) DSPC-C16FFA
python bilayersetup.py -f DSPC-34_acd16-66_1219 --DSPC 0.34 --acid16 0.66 -a 29 -r 9

# (2:1:1) DSPC-C24OH-C16FFA
python bilayersetup.py -f DSPC-50_alc24-25_acd16-25_1219 --DSPC 0.50 --alc24 0.25 --acid16 0.25 -a 32 -r 13

# (1:1) DSPC-C24OH
python bilayersetup.py -f DSPC-50_alc24-50_1219 --DSPC 0.50 --alc24 0.50 -a 32 -r 13

# (1:1) DSPC-C16FFA
python bilayersetup.py -f DSPC-50_acd16-50_1219 --DSPC 0.50 --acd16 0.50 -a 32 -r 9
