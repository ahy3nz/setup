import mbuild as mb
import os
class OH12(mb.Compound):
    def __init__(self):
        super(OH12, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH12.gro', compound=self)
