import mbuild as mb
import os
class C24OH(mb.Compound):
    def __init__(self):
        super(C24OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C24OH.gro', compound=self)
