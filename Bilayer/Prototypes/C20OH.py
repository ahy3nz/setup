import mbuild as mb
class C20OH(mb.Compound):
    def __init__(self):
        super(C20OH,self).__init__()
        mb.load('C20OH.pdb', compound=self)
        mb.translate(self, -self[0].pos)
