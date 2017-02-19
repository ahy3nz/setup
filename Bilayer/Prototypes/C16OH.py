import mbuild as mb
class C16OH(mb.Compound):
    def __init__(self):
        super(C16OH,self).__init__()
        mb.load('C16OH.pdb', compound=self)
        mb.translate(self, -self[0].pos)
