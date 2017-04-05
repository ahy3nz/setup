import mbuild as mb
import os
class FA22(mb.Compound):
    def __init__(self):
        super(FA22, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/FA22.gro', compound=self)
