import mbuild as mb
import os
class DSPC(mb.Compound):
    def __init__(self):
        super(DSPC,self).__init__()
        mb.load(os.getcwd()+'/Prototypes/new_DSPC.gro', compound=self)
        mb.translate(self, -self[0].pos)
