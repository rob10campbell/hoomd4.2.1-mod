## mpi simulation to shear a Brownian Dynamics sim of an 
## attractive colloid system 
## NOTE: requires matching GSD system file
##       for a fixed number of colloids and FIXED BOX SIZE 
## (Rob Campbell)

######### MODULE LIBRARY 
# Use HOOMD-blue
import hoomd
#import hoomd.md # molecular dynamics
# Use GSD files
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy as np
import math
import random # pseudo-random number generator
import sys
# Other
import os # miscellaneous operating system interfaces


######### CUSTOM CLASSES
#N/A


######### SIMULATION INPUTS
# General parameters
phi =  0.2 # volume fraction
rho = 3.0 # number density (per unit volume)
kT = 0.1; # system temperature
D0 = 12.0 * kT # attraction strength (gels at >=4kT)
kappa = 30.0 # range of attraction (4 (long range)- 30 (short range)), distance in BD units is approx 3/kappa

# Colloid particle details
R_C1 = 1 # 1st type colloid particle radius

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
#gamma = 6.0*np.pi*eta0*R_C1 # BD stock friction coefficient
alpha = 3.*np.pi*eta0 # BD drag coefficient, gamma=alpha*diameter in HOOMD-blue

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  print('WARNING: r_c is less than range of attraction. Increase r_c')
f_contact = 100 # magnitude of contact force (usually 100 or 1000)

# shear flow details
shear_style = 'constant'
#shear_style = 'cosinusoid' 

n_strains = 10 # number of strains or oscillation cycles
             # ex: 10 for constant strain, 5 for oscillatory, etc.
frames_per_strain = 50 # how many frames you want in each cycle
n_frames = frames_per_strain*n_strains # total simulation frames

init_velocity = False # initialize with or without a linear velocity profile

# box size
#L_X = 30 # box size in flow direction
L_Y = 30 # box size in gradient direction
#L_Z = 30 # box size in vorticity direction

# Constant
if shear_style == 'constant':
  # NOTE: reduce timestep 1e-3->1e-5 as SR increases 0.001->1.0
  dt_Integration = 1e-4 # dt! (BD timestep, may need to be smaller than DPD)

  shear_rate = 0.1
  vinf = shear_rate * L_Y # the maximum flow velocity
  tilt = 1.0 # amount of strain (desired xy tile factor 0.0 to 1.0)
  t_ramp = tilt/shear_rate # length of each strain in BD-time

else:
  # NOTE: reduce timestep 1e-3->1e-5 as SR increases 0.001->1.0
  dt_Integration = 1e-3 # dt! (BD timestep, may need to be smaller than DPD)

  Amp = 0.1  # ex: 0.005 0.01 0.05 0.1
  omega = 0.01 #2*np.pi/t_ramp # frequency (1/BD-time)
  omega_timesteps = omega*dt_Integration # frequency (1/N_timesteps)
  vinf = Amp * omega * L_Y # the maximum flow velocity
  shear_rate = Amp * omega
  t_ramp = Amp/shear_rate # length of each strain in BD-time
 
Iter_cycle = int(t_ramp / dt_Integration) # length of each shear pass (timesteps; in older code this was dt_shearing)
period= int(Iter_cycle/frames_per_strain)  # recording interval (ex: 20 frames in theta strains or omega sweeps)


# set the random seed for reproducibility
seed_value = 42


######### SIMULATION
# Checks for existing shear flow files. If none exist, begins shearing 
# from the gelation state

if os.path.exists('Shear-BD.gsd'):
  print('Shear flow file already exists. No new files created.')
  exit()
else:
  print('Shearing the Brownian Dynamics simulation with '+shear_style+' shear')

## Create a CPU simulation
device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=seed_value)

# reset timestep to zero and start simulation from the GSD file
sim.timestep=0
sim.create_state_from_gsd(filename='../3-gelation/Gelation-BD.gsd')

# get snapshot to ensure parameters only print only 
snap = sim.state.get_snapshot()

# assign particle types to groups
# (in case we want to integrate over subpopulations only,
# but would require other mods to source code)
groupB = hoomd.filter.Type(['A'])
all_colloids = hoomd.filter.Type(['A'])
all_ = hoomd.filter.Type(['A'])

#sim.state.thermalize_particle_momenta(filter=groupA, kT=kT)

# set shear flow Variant
# Constant | value
if shear_style == 'constant':

  flip = True

  shear_style_vinf = hoomd.variant.Constant(vinf)

else:

  flip = False

  # Cosinusoid | value, t_start, omega
  if shear_style == 'cosinusoid' :
    shear_style_vinf = hoomd.variant.Cosinusoid(vinf, sim.timestep, omega_timesteps)


# create neighboring list
nl = hoomd.md.nlist.Tree(buffer=0.05);

# define Morse force (attraction) interactions
morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c)

# colloid-colloid: hard particles (no deformation/overlap)
morse.params[('A', 'A')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), f_contact=f_contact)
#morse.r_cut[('B', 'B')] = (R_C1+R_C1) + (7.0/kappa) # used to assemble nl
morse.r_cut[('A','A')] = r_c+(R_C1+R_C1) # used to assemble nl


# choose integration method for the end of each timestep
# BROWNIAN (overdamped) or LANGEVIN (underdamped)
brownian=hoomd.md.methods.Brownian(filter=all_, kT=kT, alpha=alpha)
integrator=hoomd.md.Integrator(dt=dt_Integration, vinf=shear_style_vinf, forces=[morse], methods=[brownian])
#langevin = hoomd.md.methods.Langevin(filter=all_, kT=kT, alpha=alpha)
#integrator = hoomd.md.Integrator(dt=dt_Integration, vinf=shear_style_vinf, forces=[morse], methods=[langevin])
sim.operations.integrator = integrator

# set the simulation to log certain values
logger = hoomd.logging.Logger()
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
sim.operations.computes.append(thermodynamic_properties)
logger.add(thermodynamic_properties,quantities=['kinetic_temperature','pressure_tensor','virial_ind_tensor','potential_energy'])
logger.add(sim,quantities=['tps'])

# set output file
gsd_writer = hoomd.write.GSD(trigger=period, filename='Shear-BD.gsd', filter=all_, mode='wb', dynamic=['property','momentum','attribute'])
gsd_writer.write_diameter = True
#gsd_writer.maximum_write_buffer_size = 1e8 # max 100 million bytes
sim.operations.writers.append(gsd_writer)
gsd_writer.logger = logger

# shear the system 
# set the box resize style
box_resize=hoomd.update.BoxShear(trigger=1, vinf=shear_style_vinf, deltaT=dt_Integration, flip=flip)
sim.operations += box_resize

sim.run(Iter_cycle*n_strains, write_at_start=True)


if(snap.communicator.rank == 0): # if this is the start of the simulation
  print('\nNew '+shear_style+' shear Brownian Dynamics state (Shear-BD.gsd) created.')

  if shear_style == 'constant':
   print(shear_style+' shear')
   print('shear rate, SR: '+str(round(shear_rate,2)))
   print('t_ramp:', str(t_ramp))
   print('dt_Integration:', str(dt_Integration))
   print('frames per strain:', str(frames_per_strain))
   print('Iter_cycle:', str(Iter_cycle))
   print('period:', str(period))
   print('n_strains:',str(n_strains))

  else:
   print(shear_style+' shear')
   print('Amplitude: '+str(Amp,3))
   print('omega: '+str(omega,3))
   print('shear rate, SR: '+str(round(shear_rate,2)))
   print('t_ramp:', str(t_ramp))
   print('dt_Integration:', str(dt_Integration))
   print('frames per strain:', str(frames_per_strain))
   print('Iter_cycle:', str(Iter_cycle))
   print('period:', str(period))
   print('n_strains:',str(n_strains))

