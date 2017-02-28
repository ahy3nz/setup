import mbuild as mb
import os
class DSPC(mb.Compound):
    def __init__(self):
        super(DSPC, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/DSPC.gro', compound=self)


