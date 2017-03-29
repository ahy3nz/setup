import mbuild as mb
import os
class CHOL(mb.Compound):
    def __init__(self):
        super(CHOL,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/CHOL.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
