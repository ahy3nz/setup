import mbuild as mb
import os
class C22FFA(mb.Compound):
    def __init__(self):
        super(C22FFA, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C22FFA.gro', compound=self)
