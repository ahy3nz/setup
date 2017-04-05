import mbuild as mb
import os
class OH14(mb.Compound):
    def __init__(self):
        super(OH14, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH14.gro', compound=self)
