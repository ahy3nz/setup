import mbuild as mb
import os
class alc22(mb.Compound):
    def __init__(self):
        super(alc22,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/alc22.gro', compound=self)
        mb.translate(self, -self[-2].pos)
