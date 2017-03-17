import mbuild as mb
import os
class C20OH(mb.Compound):
    def __init__(self):
        super(C20OH, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C20OH.gro', compound=self)
