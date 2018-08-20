import parmed as pmd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
topfile = 'compound.top'
grofile = 'npt.gro'
temp = 305 * u.kelvin
pressure = 1 * u.bar
timestep = 2.0 * u.femtoseconds
sim_time = 100 * u.nanoseconds
n_steps = round(sim_time/timestep)

print("Reading grofiles")
top = pmd.load_file('compound.top',xyz='npt_273_400-500ns.gro')

print("Creating system from topology")
system = top.createSystem(nonbondedMethod=app.PME,
                            constraints=app.HBonds,
                            nonbondedCutoff=12.0*u.angstroms,
                            switchDistance=10.0*u.angstroms)
barostat = mm.MonteCarloMembraneBarostatp(pressure, 0.0*u.bar*u.nanometer, 
                                      temp,
                                      mm.MonteCarloMembraneBarostat.XYIsotropic,
                                      mm.MonteCarloMembraneBarostat.ZFree,
                                      100)
system.addForce(barostat)
print("Creating integrator")
integrator = mm.LangevinIntegrator(temp, 
                                    1.0/u.picoseconds, 
                                    timestep)

print("Creating Simulation")
sim = app.Simulation(top.topology, system, integrator)

print("Setting context")
sim.context.setPositions(top.positions)
sim.reporters.append(app.StateDataReporter('thermo.log', 1000, step=True, time=True,
                                            potentialEnergy=True,
                                            temperature=True,
                                            volume=True, speed=True))
sim.reporters.append(app.DCDReporter('trajectory.dcd', 5000))
sim.reporters.append(app.PDBReporter('trajectory.pdb', 5000))
sim.reporters.append(app.CheckpointReporter('trajectory.chk', 5000))
# Load the checkpoint
#with open('my_checkpoint.chk', 'rb') as f:
#        sim.context.loadCheckpoint(f.read())

print("Running MD")
sim.step(n_steps)

for reporter in sim.reporters:
    reporter.report(sim, sim.context.getState(-1))
