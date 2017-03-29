import mbuild as mb
import os
class alc12(mb.Compound):
    def __init__(self):
        super(alc12,self).__init__()
        mb.load(os.getcwd()  + '/Prototypes/alc12.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
