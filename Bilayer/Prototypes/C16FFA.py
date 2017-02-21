import mbuild as mb
class C16FFA(mb.Compound):
    def __init__(self):
        super(C16FFA,self).__init__()
        mb.load('C16FFA.pdb', compound=self)
        mb.translate(self, -self[-2].pos)
