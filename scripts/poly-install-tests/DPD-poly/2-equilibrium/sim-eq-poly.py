## mpi simulation to resolve overlaps from init.gsd
## Brings a simulation of neutral colloids suspended
## in a solvent (colloid-colloid attraction = 0)  to 
## thermal equilibrium, resolving non-physical particle
## overlaps to create a colloid suspension Equilibrium.gsd
## NOTE: for BIMODAL GELS (2 colloid types)
## NOTE: requires matching init.gsd file
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
phi = 0.20 #float(os.environ.get('phi', 0.0)) # volume fraction
percent_C1 = 1.0 #float(os.environ.get('percent_C1', 0.0)) # volume percent of type 1 colloid
percent_C2 = 0.0 #float(os.environ.get('percent_C2', 0.0)) # volume percent of type 2 colloid
poly_C1 = 0.05
poly_C2 = 0.00
rho = 3 # number density (per unit volume)

# General parameters
KT = 0.1 # system temperature
D0 = 0.0 * KT # Morse attraction strength (gels at >=4kT)
kappa = 30 # Morse attraction range (4 (long range)- 30 (short range)), distance in DPD units is approx 3/kappa 

eta0 = 0.3 # background viscosity
gamma = 4.5 # DPD controlling parameter for viscous resistance (dissipative force) 

# 375000 * 0.001 DPD times is more than enough is enough to diffuse evenly
N_time_steps = 10000 #375000 # number  of  time steps ### NOTE: maybe 200,000 timesteps is fine?
dt_Integration = 0.001 # dt! (DPD timestep)
period = 100 #2500 # recording interval

SR = 0.0 # shear parameter (shear rate = SR*L_X)

# Simulation box size (fixed by L_X)
L_X = 30 #float(os.environ.get('L_X', 0.0))
L_Y = L_X
L_Z = L_X
V_total = L_X * L_Y * L_Z # total volume of simulation box (cube)

# Solvent particle details
m_S = 1 # solvent particle mass 
R_S = 0.5 # solvent particle radius 
V_Solvents = (1 - phi) * V_total # total volume of solvent
N_Solvents = math.floor(rho * V_Solvents) # total number of solvent particles

# Colloid particle details
R_C1 = 1 #float(os.environ.get('R_C1', 0.0))  # 1st type colloid particle radius
V_C1 = (4./3.) * math.pi * R_C1 ** 3 # 1st type colloid particle volume (1 particle)
m_C1 = V_C1 * rho # 1st type colloid particle mass
V_Colloids_type1 = percent_C1 * phi * V_total # total volume of type 1 colloids
N_C1 = round(V_Colloids_type1 / V_C1) # number of 1st type of colloid particles (INT)

R_C2 = 2 #float(os.environ.get('R_C2', 0.0)) # 2nd type colloid particle radius
V_C2 = (4./3.) * math.pi * R_C2 ** 3 # 2nd type colloid particle volume (1 particle)
m_C2 = V_C2 * rho # 2nd type colloid particle mass
V_Colloids_type2 = percent_C2 * phi * V_total # total volume of type 2 colloids
N_C2 = round(V_Colloids_type2 / V_C2) # number of 2nd type of colloid particles (INT)

# colloids totals NOTE: ASSUMES 2 colloid types
N_C = N_C1 + N_C2 # total number of colloids
V_Colloids = V_Colloids_type1 + V_Colloids_type2 # total volume of all colloids

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
r0 = 0.0 # minimum inter-particle distance
f_contact = 10000.0 * KT / r_c # set colloid-colloid hard-sphere interactions 

# modified cut-off radius for colloid-solvent interactions
# NOTE: ASSUMES 2 colloid types
r_sc1_cut = (r_c**3 + R_C1**3) ** (1/3)
r_sc2_cut = (r_c**3 + R_C2**3) ** (1/3)

# Total number of particles in the simulation
N_total = int(N_Solvents + N_C)


######### SIMULATION
## Checks for existing equilibrium files. If none are found, brings the 
## initial random distribution of particles to thermal equilibrium. 

if os.path.exists('Equilibrium-poly-DPD.gsd'):
  print("Equilibrium file already exists. No new files created.")
else:
  print("Initialization state is being brought to equilibrium")
  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations

  # start the simulation from the initialized system
  sim.timestep = 0 # set initial timestep to 0
  sim.create_state_from_gsd(filename='../1-initialize/init-poly-DPD.gsd')
  # optional: force domain decomposition to predetermined viable size 3x3x3=27cores
  #sim.create_state_from_gsd(filename='../1-initialize/init.gsd', domain_decomposition=(5,5,5))


  # assign particle types to groups 
  # (in case we want to integrate over subpopulations only, 
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['A'])
  groupB = hoomd.filter.Type(['B'])
  groupC = hoomd.filter.Type(['C'])
  all_ = hoomd.filter.Type(['A','B','C'])
  all_colloids = hoomd.filter.Type(['B','C'])

  # thermalize (aka bring to thermal equilibrium) the system
  sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define DPD Morse force (attraction) interactions
  morse = hoomd.md.pair.DPDMorse(nlist=nl, kT=KT, default_r_cut=1.0 * r_c, bond_calc=False, poly=poly_C1+poly_C2)

  # solvent-solvent: soft particles (allow deformation/overlap)
  morse.params[('A','A')] = dict(A0=25.0 * KT / r_c, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=0.0, rcut=r_c, poly=0.0) # force calc
  morse.r_cut[('A','A')] = r_c # used to assemble nl

  # solvent-colloid: soft solvent particles (allow deformation/overlap)
  morse.params[('A','B')] = dict(A0=25.0 * KT / r_sc1_cut, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=R_C1, rcut=r_sc1_cut - (0.0 + R_C1), poly=0.0+poly_C1) # force calc
  morse.r_cut[('A','B')] = r_sc1_cut # used to assemble nl

  # solvent-colloid: soft solvent particles (allow deformation/overlap)
  morse.params[('A','C')] = dict(A0=25.0 * KT / r_sc2_cut, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=R_C2, rcut=r_sc2_cut - (0 + R_C2), poly=0.0+poly_C2) # force calc
  morse.r_cut[('A','C')] = r_sc2_cut # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','B')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C1, a2=R_C1, rcut=r_c, poly=poly_C1) # force calc
  morse.r_cut[('B','B')] = (r_c + R_C1 + R_C1) # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','C')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C1, a2=R_C2, rcut=r_c, poly=poly_C1+poly_C2) # force calc
  morse.r_cut[('B','C')] = (r_c + R_C1 + R_C2) # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('C','C')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C2, a2=R_C2, rcut=r_c, poly=poly_C2) # force calc
  morse.r_cut[('C','C')] = (r_c + R_C2 + R_C2) # used to assemble nl

  # choose integration method for the end of each timestep
  nve = hoomd.md.methods.ConstantVolume(filter=all_, thermostat=None)
  integrator=hoomd.md.Integrator(dt=dt_Integration, SR=SR*L_Y, forces=[morse], methods=[nve])
  sim.operations.integrator = integrator

  # set the simulation to log certain values
  logger = hoomd.logging.Logger()
  thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
  sim.operations.computes.append(thermodynamic_properties)
  logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 
    'pressure_tensor', 'virial_ind_tensor', 'potential_energy'])
  logger.add(sim, quantities=['tps'])
  logger.add(morse, quantities=['forces','virials'])

  # save outputs
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Equilibrium.gsd", 
    filter=all_, mode='wb', dynamic=['property','momentum','attribute'])
  gsd_writer.write_diameter = True
  sim.operations.writers.append(gsd_writer)
  gsd_writer.logger = logger

  '''
  # save colloids ONLY for smaller file-size
  # NOTE: you MUST save all particles to use the output to run another sim (gel, shear, etc.)
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Equilibrium_Colloids.gsd", 
    filter=all_colloids, mode='wb', dynamic=['property','momentum','attribute'])
  sim.operations.writers.append(gsd_writer)
  gsd_writer.log = logger
  '''

  # run simulation!
  sim.run(N_time_steps, write_at_start=True)

  print("New equilibrium state (Equilibrium-poly-DPD.gsd) created.")
