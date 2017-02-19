import mbuild as mb
class C24OH(mb.Compound):
    def __init__(self):
        super(C24OH,self).__init__()
        mb.load('C24OH.pdb', compound=self)
        mb.translate(self, -self[0].pos)
