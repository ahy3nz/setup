import mbuild as mb
class CHOL(mb.Compound):
    def __init__(self):
        super(CHOL,self).__init__()
        mb.load('CHOL.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
