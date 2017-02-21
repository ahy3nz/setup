import mbuild as mb
class ISIS(mb.Compound):
    def __init__(self):
        super(ISIS,self).__init__()
        mb.load('ISIS.pdb', compound=self)
        mb.translate(self, -self[0].pos)
