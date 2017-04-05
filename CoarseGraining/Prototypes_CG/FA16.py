import mbuild as mb
import os
class FA16(mb.Compound):
    def __init__(self):
        super(FA16, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/FA16.pdb', compound=self)
