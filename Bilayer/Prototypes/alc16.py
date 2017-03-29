import mbuild as mb
import os
class alc16(mb.Compound):
    def __init__(self):
        super(alc16,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc16.gro', compound=self)
        mb.translate(self, -self[-2].pos)
