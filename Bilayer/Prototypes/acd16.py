import mbuild as mb
import os
class acd16(mb.Compound):
    def __init__(self):
        super(acd16,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/acd16.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
