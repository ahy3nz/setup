import mbuild as mb
import os
class w(mb.Compound):
    def __init__(self):
        super(w, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/w.gro', compound=self)

