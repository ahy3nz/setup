import mbuild as mb
import os
class OH22(mb.Compound):
    def __init__(self):
        super(OH22, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH22.pdb', compound=self)
