## mpi simulation to bring a Brownian Dynamics sim of an 
## attractive colloid system from an equilibrium suspension 
## to a quasi-steady state gel
## NOTE: requires matching Equilibrium-BD.gsd file
##       for a FIXED NUMBER OF COLLOIDS and variable box size
## (Rob Campbell)


######### MODULE LIBRARY 
# Use HOOMD-blue
import hoomd
import hoomd.md # molecular dynamics
# Use GSD files
import gsd # provides Python API; MUST import sub packages explicitly (as below)
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy
import math
import random # pseudo-random number generator
# Other
import os # miscellaneous operating system interfaces


######### CUSTOM CLASSES
#N/A


######### SIMULATION INPUTS
# General parameters
phi =  0.2 # volume fraction
rho = 3.0 # number density (per unit volume)
KT = 0.1; # system temperature
D0 = 12.0 * KT # attraction strength (gels at >=4kT)

N_time_steps = 1500000 # number of time steps
dt_Integration = 0.001 # dt! (BD timestep, may need to be smaller than DPD)
period = 10000 # recording interval

# Colloid particle details
N_C1 = 1000 # number of 1st type of colloidal particles
R_C1 = 1  # 1st type colloid particle radius
V_C1 = (4./3.) * math.pi * R_C1 ** 3 # 1st type colloid particle volume (1 particle)
m_C1 = V_C1 * rho # 1st type colloid particle mass

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
gamma = 6.0*numpy.pi*eta0*R_C1 # BD stock friction coefficient

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  print('WARNING: r_c is less than range of attraction. Increase r_c')
r0 = 0.0 # minimum inter-particle distance
kappa = 30.0 # range of attraction (4 (long range)- 30 (short range)), distance in BD units is approx 3/kappa
f_contact = 10000.0 * KT / r_c # set colloid-colloid hard-sphere interactions 
#bond_calc = False # do you want to track what bonds form and break? True=yes, False=no

# Simulation box size (NOTE: calculated from # colloids)
L_X = (N_C1 * V_C1 / phi)**(1./3.)
L_Y = L_X
L_Z = L_X 

# Volumes
V_total = L_X*L_Y*L_Z; # total volume of simulation box (cube)
# NOTE: volume calculation currently ASSUMES 1 colloid type and 1 solvent type
V_Colloids = N_C1 * V_C1; # total volume of all colloid particles

# Total number of particles in the simulation
N_total = int(N_C1)


######### SIMULATION
# Checks for existing gelation files. If none exist, continues from
# the equilibrium state towards a quasi-steady state gel (no shear).

if os.path.exists('Gelation-BD.gsd'):
  print('Gelation file already exists. No new files created.')
else:
  print('Brownian Dynamics equilibrium state is being brought to quasi-steady state gelation')
  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations

  # start the simulation from Equilibrium (and don't set timestep to 0)
  sim.create_state_from_gsd(filename='../2-equilibrium/Equilibrium-BD.gsd')

  # assign particle types to groups
  # (in case we want to integrate over subpopulations only,
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['A'])
  all_ = hoomd.filter.Type(['A'])

  # don't want to thermalize the system after equilibrium because
  # we are tracking velocity during gelation
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define Morse force (attraction) interactions
  # TO BE ADDED: bond_calc=bond_calc
  morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c)

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('A','A')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1))	
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
  logger.add(thermodynamic_properties,quantities=['kinetic_temperature',
		'pressure_tensor','virial_ind_tensor','potential_energy'])
  logger.add(sim, quantities=['tps'])

  # set output file
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Gelation-BD.gsd", 
    filter=all_, mode='wb', dynamic=['property','momentum','attribute'])

  # don't need to thermalize the system after equilibrium because
  # we are tracking velocity during gelation
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)
  gsd_writer.maximum_write_buffer_size = period

  # save outputs
  sim.operations.writers.append(gsd_writer)
  gsd_writer.logger = logger

  # run simulation!
  # (and write the initial state (e.g. the last frame of Equilibrium) in this file!)
  sim.run(N_time_steps, write_at_start=True)
	
  print('New Brownian Dynamics gelation state (Gelation-BD.gsd) created.')
