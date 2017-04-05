import mbuild as mb
import os
class OH20(mb.Compound):
    def __init__(self):
        super(OH20, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH20.pdb', compound=self)
