import mbuild as mb
import os
class acd12(mb.Compound):
    def __init__(self):
        super(acd12,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/acd12.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
