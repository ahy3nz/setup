import mbuild as mb
import os
class C12OH(mb.Compound):
    def __init__(self):
        super(C12OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C12OH.gro', compound=self)
