import mbuild as mb
import os
class C14OH(mb.Compound):
    def __init__(self):
        super(C14OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C14OH.gro', compound=self)
