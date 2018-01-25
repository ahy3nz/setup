import mbuild as mb
import os
class PCN(mb.Compound):
    def __init__(self):
        super(PCN,self).__init__()
        mb.load(os.getcwd()+'/Prototypes/PCN.gro', compound=self)
        #mb.translate(self, -self[0].pos)
        self.translate(-self[0].pos)
