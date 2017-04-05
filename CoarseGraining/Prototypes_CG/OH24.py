import mbuild as mb
import os
class OH24(mb.Compound):
    def __init__(self):
        super(OH24, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH24.pdb', compound=self)
