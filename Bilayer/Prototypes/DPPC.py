import mbuild as mb
import os
class DPPC(mb.Compound):
    def __init__(self):
        super(DPPC,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/DPPC.pdb', compound=self)
        mb.translate(self, -self[0].pos)
