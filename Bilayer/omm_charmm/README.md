# Running CHARMM in OpenMM
* What's mainly important is the longrange treatment and cutoffs. 
* There are lots of ways to create the openMM `System` object (Foyer, Parmed). 
Using Parmed is easiest because it can read in gromacs files to generate the appropriate
`System` object
* OpenMM has `StateDataReporter`, `PDBReporter`, `CheckpointReporter`, `DCDReporter` -
use them.
* Note the use of the `MonteCarloMembraneBarostat`
* Realize that `PDBReporter` might not carry over box vectors. This information is
stored in the `sim.context.state`, however.
