import mbuild as mb
import os
class H2O(mb.Compound):
    def __init__(self):
        super(H2O,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/H2O.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
