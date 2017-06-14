import mbuild as mb
import os
class HOH(mb.Compound):
    def __init__(self):
        super(HOH,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/HOH.pdb', compound=self, use_parmed=True)
        mb.translate(self, -self[-1].pos)
