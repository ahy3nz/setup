import mbuild as mb
import os
class alc18(mb.Compound):
    def __init__(self):
        super(alc18,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc18.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
