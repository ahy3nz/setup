import mbuild as mb
import os
class water(mb.Compound):
    def __init__(self):
        super(water, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/water.gro', compound=self)

