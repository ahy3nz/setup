import mbuild as mb
import os
class alc24(mb.Compound):
    def __init__(self):
        super(alc24,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc24.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
