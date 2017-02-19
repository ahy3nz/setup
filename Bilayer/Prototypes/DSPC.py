import mbuild as mb
class DSPC(mb.Compound):
    def __init__(self):
        super(DSPC,self).__init__()
        mb.load('DSPC.pdb', compound=self)
        mb.translate(self, -self[0].pos)
