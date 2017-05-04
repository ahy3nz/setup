import mbuild as mb
import os
class acd20(mb.Compound):
    def __init__(self):
        super(acd20,self).__init__()
        mb.load(os.getcwd()  + '/Prototypes/acd20.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
