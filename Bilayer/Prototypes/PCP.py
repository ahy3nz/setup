import mbuild as mb
import os
class PCP(mb.Compound):
    def __init__(self):
        super(PCP,self).__init__()
        mb.load(os.getcwd()+'/Prototypes/PCP.gro', compound=self)
        #mb.translate(self, -self[0].pos)
        self.translate(-self[0].pos)
