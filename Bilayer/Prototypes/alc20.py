import mbuild as mb
import os
class alc20(mb.Compound):
    def __init__(self):
        super(alc20,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc20.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
