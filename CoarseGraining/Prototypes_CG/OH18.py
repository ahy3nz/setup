import mbuild as mb
import os
class OH18(mb.Compound):
    def __init__(self):
        super(OH18, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/OH18.gro', compound=self)
