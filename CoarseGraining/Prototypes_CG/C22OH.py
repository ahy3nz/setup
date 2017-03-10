import mbuild as mb
import os
class C22OH(mb.Compound):
    def __init__(self):
        super(C22OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C22OH.gro', compound=self)
