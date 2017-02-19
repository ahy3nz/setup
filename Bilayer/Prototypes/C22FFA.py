import mbuild as mb
class C22FFA(mb.Compound):
    def __init__(self):
        super(C22FFA,self).__init__()
        mb.load('C22FFA.pdb', compound=self)
        mb.translate(self, -self[0].pos)
