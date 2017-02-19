import mbuild as mb
class C18OH(mb.Compound):
    def __init__(self):
        super(C18OH,self).__init__()
        mb.load('C18OH.pdb', compound=self)
        mb.translate(self, -self[0].pos)
