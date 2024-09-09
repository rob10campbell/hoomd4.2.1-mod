## mpi simulation to bring a Brownian Dynamics sim of an 
## attractive colloid system (with added angular repulsion in 
## many-body interactions of certain angles) from an 
## equilibrium suspension to a quasi-steady state gel
## NOTE: requires matching Equilibrium-BD.gsd file
##       for a fixed number of colloids and FIXED BOX SIZE 
## (Rob Campbell)


######### MODULE LIBRARY 
# Use HOOMD-blue
import hoomd
import hoomd.md # molecular dynamics
# Use GSD files
import gsd # provides Python API; MUST import sub packages explicitly (as below)
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy as np
import math
import random # pseudo-random number generator
# Other
import os # miscellaneous operating system interfaces

#pip install mpi4py
from mpi4py import MPI # for tracking MPI processes
# track MPI rank and comm.size (1 if serial)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if comm.size > 1:
  print("MPI rank: "+str(rank))


######### CUSTOM CLASSES
#N/A


######### SIMULATION INPUTS
# General parameters
phi =  0.1 # volume fraction
rho = 3.0 # number density (per unit volume)
KT = 0.1 # system temperature

# Morse
D0 = 12.0 * KT # attraction strength (gels at >=4kT)
kappa = 30.0 # range of attraction (4 (long range)- 30 (short range)), distance in BD units is approx 3/kappa

# angular repulsion params
theta_bar = 65 # angle in degrees
theta_bar_rad = theta_bar * np.pi / 180 # angle in radians
B = 0.0 #1 * KT #67.27 # magnitude of angular repulsion
w = 0.30 # width of the angular repulsion potential

N_time_steps = 10000 #1500000 # number of time steps
dt_Integration = 0.001 # dt! (BD timestep, may need to be smaller than DPD)
period = 100 #10000 # recording interval

# Colloid particle details
R_C1 = 1 # 1st type colloid particle radius

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
#gamma = 6.0*np.pi*eta0*R_C1 # BD stock friction coefficient
alpha = 3.*np.pi*eta0 # BD drag coefficient, gamma=alpha*diameter in HOOMD-blue

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  if rank == 0 or comm.size == 1:
    print('WARNING: r_c is less than range of attraction. Increase r_c')
r0 = 0.0 # minimum inter-particle distance
f_contact = 0 #100 # magnitude of contact force (usually 100 or 1000)
#bond_calc = False # do you want to track what bonds form and break? True=yes, False=no

seed_value = 45

######### SIMULATION
# Checks for existing gelation files. If none exist, continues from
# the equilibrium state towards a quasi-steady state gel (no shear).

if os.path.exists('Gelation-BD.gsd'):
  if rank == 0 or comm.size == 1:
    print('Gelation file already exists. No new files created.')
else:
  if rank == 0 or comm.size == 1:
    print('Brownian Dynamics equilibrium state is being brought to quasi-steady state gelation')

  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=seed_value) # set seed to a fixed value for reproducible simulations

  # start the simulation from Equilibrium (and don't set timestep to 0)
  sim.create_state_from_gsd(filename='../../bond-limits/L30_phi10/2-equilibrium/Equilibrium-BD.gsd')

  # assign particle types to groups
  # (in case we want to integrate over subpopulations only,
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['B'])
  all_ = hoomd.filter.Type(['B'])

  # don't want to thermalize the system after equilibrium because
  # we are tracking velocity during gelation
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define Morse force (attraction) interactions
  # TO BE ADDED: bond_calc=bond_calc
  morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c)
  #morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c, B=B, theta_bar=theta_bar_rad)
  #angle = hoomd.md.many_body.AngularRepulsion(default_r_cut=3.0, nlist=nl)

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','B')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), f_contact=f_contact)
  morse.r_cut[('B','B')] = r_c+(R_C1+R_C1) # used to assemble nl

  #angle.params[('B','B')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), f_contact=f_contact)
  #angle.params[('B','B')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), theta_bar=theta_bar_rad, B=B, w=w, scaled_D0=False)
  #angle.r_cut[('B','B')] = r_c+(R_C1+R_C1) # used to assemble nl

  #tersoff.params[('B', 'B')] = dict(magnitudes=(2.0, 1.0), lambda3=5.0)

  # choose integration method for the end of each timestep
  # BROWNIAN (overdamped) or LANGEVIN (underdamped)
  brownian = hoomd.md.methods.Brownian(filter=all_, kT=KT, alpha=alpha)
  #integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[morse,tersoff], methods=[brownian])
  integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[morse], methods=[brownian])
  #integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[angle], methods=[brownian])
  #langevin = hoomd.md.methods.Langevin(filter=all_, kT=KT, alpha=alpha)
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
	
  if rank == 0 or comm.size == 1:
    print('New Brownian Dyanmics gelation state (Gelation-BD.gsd) created.')
