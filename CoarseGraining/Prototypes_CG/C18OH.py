import mbuild as mb
import os
class C18OH(mb.Compound):
    def __init__(self):
        super(C18OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C18OH.gro', compound=self)
