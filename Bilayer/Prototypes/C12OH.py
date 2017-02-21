import mbuild as mb
class C12OH(mb.Compound):
    def __init__(self):
        super(C12OH,self).__init__()
        mb.load('C12OH.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
