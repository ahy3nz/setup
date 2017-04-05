import mbuild as mb
import os
class OH16(mb.Compound):
    def __init__(self):
        super(OH16, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH16.gro', compound=self)
