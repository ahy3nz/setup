import mbuild as mb
import os
class ISIS(mb.Compound):
    def __init__(self):
        super(ISIS,self).__init__()
        mb.load(os.getcwd() + '/Prototypes/new_ISIS.gro', compound=self)
        mb.translate(self, -self[0].pos)
