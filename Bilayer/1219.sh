#!/bin/bash
# (1:1) DSPC-C12OH
python bilayersetup.py -f DSPC-50_alc12-50_1219 --DSPC 0.50 --alc22 0.50 -a 31 -r 9

# (1:4) DSPC-C24OH
python bilayersetup.py -f DSPC-20_alc24-80_1219 --DSPC 0.20 --alc24 0.80 -a 26 -r 20

# (1:2) DSPC-C12OH

# (1:2) DSPC-C22FFA
python bilayersetup.py -f DSPC-34_acd22-66_1219 --DSPC 0.34 --acd22 0.66 -a 29 -r 18

# (1:1:1) DSPC-C22OH-C22FFA
python bilayersetup.py -f DSPC-34_alc22-33_acd22-33_1219 --DSPC 0.34 --alc22 0.33 --acd22 0.33 -a 24 -r 6

# (1:1:1) DSPC-C24OH-C16FFA

# (1:2) DSPC-C24OH

# (1:2) DSPC-C16FFA

# (2:1:1) DSPC-C24OH-C16FFA

# (1:1) DSPC-C24OH

# (1:1) DSPC-C16FFA
