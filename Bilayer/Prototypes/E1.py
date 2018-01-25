import mbuild as mb
import os
class E1(mb.Compound):
    def __init__(self):
        super(E1,self).__init__()
        mb.load(os.getcwd()+'/Prototypes/E1.gro', compound=self)
        #mb.translate(self, -self[0].pos)
        self.translate(-self[0].pos)
