import mbuild as mb
import os
class alc14(mb.Compound):
    def __init__(self):
        super(alc14,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc14.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
