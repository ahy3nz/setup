import mbuild as mb
import os
class PC(mb.Compound):
    def __init__(self):
        super(PC,self).__init__()
        mb.load(os.getcwd()+'/Prototypes/PC.gro', compound=self)
        #mb.translate(self, -self[0].pos)
        self.translate(-self[0].pos)
