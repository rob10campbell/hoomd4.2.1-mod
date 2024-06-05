## mpi simulation to bring a DPD sim of an attractive 
## colloid system from equilibrium suspension to a 
## quasi-steady state gel
## NOTE: for variable colloid size (POLYDISPERSITY)
## NOTE: requires matching Equilibrium-DPD.gsd file
##       for a fixed number of colloids and FIXED BOX SIZE 
## (Rob Campbell)


#########  MODULE LIBRARY
# Use HOOMD-blue
import hoomd
import hoomd.md # molecular dynamics
# Use GSD files
import gsd # provides Python API; MUST import sub packages explicitly (as below)
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy
import math
import random # psuedo-random number generator
# Other
import os # miscellaneous operating system interfaces


######### CUSTOM CLASSES
#N/A


######### SIMULATION INPUTS
# Density parameters
phi =  0.2 # volume fraction
poly = 0.05  # standard deviation in size (AKA % polydispersity)
rho = 3 # number density (per unit volume)

# General parameters
KT = 0.1 # system temperature
D0 = 12.0 * KT # attraction strength (gels at >=4kT)
scaled_D0 = False
kappa = 30 # range of attraction (4 (long range)- 30 (short range)), distance in DPD units is approx 3/kappa

eta0 = 0.3 # background viscosity
gamma = 4.5 # DPD controlling parameter for viscous resistance (dissipative force) 

N_time_steps = 2000000 # number  of  time steps
dt_Integration = 0.001 # dt! (DPD timestep)
period = 10000 # recording interval

# Colloid particle details
R_C1 = 1 # 1st type colloid particle radius

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter (r_cut = # * r_c) 
r0 = 0.0 # minimum inter-particle distance
f_contact = 10000.0 * KT / r_c # set colloid-colloid hard-sphere interactions 

# modified cut-off radius for colloid-solvent interactions
# NOTE: ASSUMES 2 colloid types
r_sc1_cut = (r_c**3 + R_C1**3) ** (1/3)


######### SIMULATION
## Checks for existing equilibrium files. If none are found, brings the 
## initial random distribution of particles to thermal equilibrium. 

if os.path.exists('Gelation-yly-DPD.gsd'):
  print("Gelation file already exists. No new files created.")
else:
  print("Polydisperse DPD equilibrium state is continuing to be brought to quasi-steady state gelation")
  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations

  # start the simulation from Equilibrium (and don't reset timestep to 0)
  sim.create_state_from_gsd(filename='../2-equilibrium/Equilibrium-poly-DPD.gsd')
  # optional: force domain decomposition to predetermined viable size 3x3x3=27cores
  #sim.create_state_from_gsd(filename='../../../2-equilibrium/Equilibrium-poly-DPD.gsd', domain_decomposition=(3,3,3))

  # assign particle types to groups 
  # (in case we want to integrate over subpopulations only, 
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['A'])
  groupB = hoomd.filter.Type(['B'])
  all_ = hoomd.filter.Type(['A','B'])
  all_colloids = hoomd.filter.Type(['B'])

  # don't need to thermalize the system after equilibrium because
  # we are tracking velocity during gelation
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define DPD Morse force (attraction) interactions
  morse = hoomd.md.pair.DPDMorse(nlist=nl, kT=KT, default_r_cut=1.0 * r_c, bond_calc=False)

  # solvent-solvent: soft particles (allow deformation/overlap)
  morse.params[('A','A')] = dict(A0=25.0 * KT / r_c, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=0.0, rcut=r_c) # force calc
  morse.r_cut[('A','A')] = r_c # used to assemble nl

  # solvent-colloid: soft solvent particles (allow deformation/overlap)
  morse.params[('A','B')] = dict(A0=25.0 * KT / r_sc1_cut, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=R_C1, rcut=r_sc1_cut - (0.0 + R_C1)) # force calc
  morse.r_cut[('A','B')] = r_sc1_cut # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','B')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C1, a2=R_C1, rcut=r_c, scaled_D0=scaled_D0) # force calc
  morse.r_cut[('B','B')] = (r_c + R_C1 + R_C1) # used to assemble nl

  # choose integration method for the end of each timestep
  nve = hoomd.md.methods.ConstantVolume(filter=all_, thermostat=None)
  integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[morse], methods=[nve])
  sim.operations.integrator = integrator

  # set the simulation to log certain values
  logger = hoomd.logging.Logger()
  thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
  sim.operations.computes.append(thermodynamic_properties)
  logger.add(thermodynamic_properties,quantities=['kinetic_temperature',
    'pressure_tensor','virial_ind_tensor','potential_energy'])
  logger.add(sim, quantities=['tps'])

  # set output file
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Gelation-poly-DPD.gsd", 
    filter=all_, mode='wb', dynamic=['property','momentum','attribute'])

  # save diameters
  gsd_writer.write_diameter = True

  # [optional] set buffer size (how often data is saved to the file, large buffers can increase performace, but can lead to lost data if sim is cancelled or times-out)
  #gsd_writer.maximum_write_buffer_size = 1e8 # max 100 million bytes

  # save outputs
  sim.operations.writers.append(gsd_writer)
  gsd_writer.logger = logger

  # save colloids ONLY for smaller file-size
  # NOTE: you MUST save all particles to use the output to run another sim (gel, shear, etc.)
  gsd_writer_colloids = hoomd.write.GSD(trigger=period, filename="Gelation_Colloids-poly-DPD.gsd", 
    filter=all_colloids, mode='wb', dynamic=['property','momentum','attribute'])
  gsd_writer_colloids.write_diameter = True
  #gsd_writer_colloids.maximum_write_buffer_size = 1e8 # max 100 million bytes
  sim.operations.writers.append(gsd_writer_colloids)
  gsd_writer_colloids.logger = logger
  
  # run simulation!
  # (and write the initial state (e.g. the last frame of Equilibrium) in this file!)
  sim.run(N_time_steps, write_at_start=True)
	
  print('New polydipserse DPD gelation state (Gelation-poly-DPD.gsd) created.')
