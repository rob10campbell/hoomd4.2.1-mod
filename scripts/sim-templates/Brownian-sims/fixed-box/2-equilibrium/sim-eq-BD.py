## mpi simulation to resolve overlaps from init-BD.gsd
## Brings a Brownian Dynamics simulation of attractive colloids
## to thermal equilibrium and creates Equilibrium-BD.gsd
## NOTE: requires matching init-BD.gsd file
##       for a variable number of colloids and FIXED BOX SIZE 
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
# General parameters
phi =  0.2 # volume fraction
rho = 3 # number density (per unit volume)
KT = 0.1 # system temperature
D0 = 0.0 * KT # attraction strength (gels at >=4kT)
kappa = 30.0 # range of attraction (4 (long range)- 30 (short range)), distance in BD units is approx 3/kappa 

N_time_steps = 100000 # number  of  time steps
dt_Integration = 0.001 # dt! (BD timestep, may need to be smaller than DPD)
period = 10000 # recording interval

# Colloid particle details
R_C1 = 1 # 1st type colloid particle radius

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
gamma = 6.0*numpy.pi*eta0*R_C1 # BD stock friction coefficient

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  print('WARNING: r_c is less than range of attraction. Increase r_c')
f_contact = 100 # magnitude of contact force (usually 100 or 1000)



######### SIMULATION
## Checks for existing equilibrium files. If none are found, brings the 
## initial random distribution of particles to thermal equilibrium. 

if os.path.exists('Equilibrium-BD.gsd'):
  print("Equilibrium file already exists. No new files created.")
else:
  print("Brownian Dynamics initialization state is being brought to equilibrium")
  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations

  # start the simulation from the initialized system
  sim.timestep = 0 # set initial timestep to 0
  sim.create_state_from_gsd(filename='../1-initialize/init-BD.gsd')

  # assign particle types to groups 
  # (in case we want to integrate over subpopulations only, 
  # but would require other mods to source code)
  groupB = hoomd.filter.Type(['A'])
  all_ = hoomd.filter.Type(['A'])

  # thermalize (aka bring to thermal equilibrium) the system
  sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define Morse force (attraction) interactions
  morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c)

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('A','A')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), f_contact=f_contact)	
  morse.r_cut[('A','A')] = r_c+(R_C1+R_C1) # used to assemble nl

  # choose integration method for the end of each timestep
  # BROWNIAN (overdamped) or LANGEVIN (underdamped)
  brownian = hoomd.md.methods.Brownian(filter=all_, kT=KT, default_gamma=gamma)
  integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[morse], methods=[brownian])
  #langevin = hoomd.md.methods.Langevin(filter=all_, kT=KT, default_gamma=gamma)
  #integrator = hoomd.md.Integrator(dt=dt_Integration, forces=[morse], methods=[langevin])
  sim.operations.integrator = integrator

  # set the simulation to log certain values
  logger = hoomd.logging.Logger()
  thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
  sim.operations.computes.append(thermodynamic_properties)
  logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 
    'pressure_tensor', 'virial_ind_tensor', 'potential_energy'])
  logger.add(sim, quantities=['tps'])

  # set output file
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Equilibrium-BD.gsd", 
    filter=all_, mode='wb', dynamic=['property','momentum','attribute'])

  # save diameters
  gsd_writer.write_diameter = True

  # [optional] set buffer size (how often data is saved to the file, large buffers can increase performace, but can lead to lost data if sim is cancelled or times-out)
  #gsd_writer.maximum_write_buffer_size = 1e8 # max 100 million bytes

  # save outputs
  sim.operations.writers.append(gsd_writer)
  gsd_writer.logger = logger

  # run simulation!
  # (and write the initial state (e.g. the last frame of Equilibrium) in this file!)
  sim.run(N_time_steps, write_at_start=True)

  print("New Brownian Dynamics equilibrium state (Equilibrium-BD.gsd) created.")
