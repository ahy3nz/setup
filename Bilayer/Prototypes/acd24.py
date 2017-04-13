import mbuild as mb
import os
class acd24(mb.Compound):
    def __init__(self):
        super(acd24,self).__init__()
        mb.load(os.getcwd()  + '/Prototypes/acd24.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
