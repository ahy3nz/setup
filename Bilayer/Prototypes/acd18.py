import mbuild as mb
import os
class acd18(mb.Compound):
    def __init__(self):
        super(acd18,self).__init__()
        mb.load(os.getcwd()  + '/Prototypes/acd18.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
