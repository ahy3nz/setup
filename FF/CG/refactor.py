import numpy as np
import pdb
import os
import subprocess

refactor_dict = {"Q0": "PCN", "Qa": "PCP",
        "Na": "E1", "C1":"C3", "C2":"C2",
        "P3":"P3",
        "P4":"P4",
        "Nda":"Nda"}

all_files = [thing for thing in os.listdir() if '.txt' in thing and os.path.isfile(thing)]
for textfile in all_files:
    old_name = textfile
    new_name = textfile
    for key, val in refactor_dict.items():
        new_name = new_name.replace(key, val)
    print("replacing {} with {}".format(old_name, new_name))
    p=subprocess.Popen('cp {} {}'.format(old_name,new_name), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
