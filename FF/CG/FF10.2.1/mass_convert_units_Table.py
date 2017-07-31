import os
import subprocess
import pdb

with open('out.log','w') as outfile: 
    for f in os.listdir():
        if '.txt' in f:
            print(f)
            p = subprocess.Popen('python convert_units_Table.py -f {}'.format(f), shell=True, stdout=outfile, stderr=outfile)
            p.wait()
