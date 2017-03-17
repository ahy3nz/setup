import mbuild as mb
import os
class C16FFA(mb.Compound):
    def __init__(self):
        super(C16FFA, self).__init__()
        mb.load(os.getcwd()+'/Prototypes_CG/C16FFA.gro', compound=self)
