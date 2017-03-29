import mbuild as mb
import os
class acd22(mb.Compound):
    def __init__(self):
        super(acd22,self).__init__()
        mb.load(os.getcwd()  + '/Prototypes/acd22.gro', compound=self)
        mb.translate(self, -self[-2].pos)
