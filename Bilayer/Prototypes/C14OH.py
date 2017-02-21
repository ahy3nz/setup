import mbuild as mb
class C14OH(mb.Compound):
    def __init__(self):
        super(C14OH,self).__init__()
        mb.load('C14OH.pdb', compound=self)
        mb.translate(self, -self[-1].pos)
