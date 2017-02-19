import mbuild as mb
class DPPC(mb.Compound):
    def __init__(self):
        super(DPPC,self).__init__()
        mb.load('DPPC.pdb', compound=self)
        mb.translate(self, -self[0].pos)
