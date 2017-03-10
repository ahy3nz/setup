import mbuild as mb
import os
class C16OH(mb.Compound):
    def __init__(self):
        super(C16OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C16OH.gro', compound=self)
