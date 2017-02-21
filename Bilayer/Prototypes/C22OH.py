import mbuild as mb
class C22OH(mb.Compound):
    def __init__(self):
        super(C22OH,self).__init__()
        mb.load('C22OH.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
