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

#N_time_steps = 15000000 # number of time steps
dt_Integration = 0.0001 # dt! (BD timestep, may need to be smaller than DPD)
period = 100000 # recording interval

# Colloid particle details
R_C1 = 1 # 1st type colloid particle radius

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
gamma = 6.0*numpy.pi*eta0*R_C1 # BD stock friction coefficient

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if rank == 0 or comm.size == 1:
  if r_c < (3/kappa):
    print('WARNING: r_c is less than range of attraction. Increase r_c')
r0 = 0.0 # minimum inter-particle distance
f_contact = 100 # magnitude of contact force (usually 100 or 1000)
#bond_calc = False # do you want to track what bonds form and break? True=yes, False=no


# shear flow details
# set shear_style
shear_style = 'constant' 
shear_style = 'oscillatory'
shear_style = 'sinusoid' 
shear_style = 'cosinusoid' 

t_ramp = int(1e5) # length of each shear pass 

#L_X = 30 # box size in flow direction
L_Y = 30 # box size in gradient direction
#L_Z = 30 # box size in vorticity direction

# Constant
if shear_style == 'constant':

  shear_rate = 0.01
  vinf = shear_rate * L_Y

  """
  # What are the times?
  """
  # old code:
  theta = 1.0 # xy tilt factor
  delta_T_shearing = round(theta/SR/dt_Integration) # timestep for shear
  nframe_strain = delta_T_shearing/10 # 20 frames in theta strains
  # --> run n_strains times for delta_T_shearing+1 timesteps

  # t_ramp = 1e5 = dT_shearing/10 = nframe_strain = period?
  dT_shearing = t_ramp*10 # length of each shear pass (ex: 20 frames in theta strains)
  theta = dT_shearing*dt_Integration*shear_rate # desired xy tilt factor (0.0 to 1.0)

  n_cycles = 10 # number of strains

  ###shear_style_vinf = hoomd.variant.Constant(vinf)
  ###shear_style_erate = hoomd.variant.Cosinusoid(shear_rate, sim.timestep, dt_Integration)

  # t_ramp is the period? it must be the number of timesteps per interval
  # t_ramp / dt_Integration = ntimesteps for a ramp
  # save every 100 timesteps?
  Iter_cycle = int(t_ramp / dt_Integration) # DPD times
  period = int(Iter_cycle/100) #

else:
  Amp = 0.1 
  omega = 2*np.pi/t_ramp # frequency (1/recording-interval)
  vinf = Amp * omega * L_Y
  shear_rate = Amp * omega

  Iter_cycle = int(2*np.pi/omega / dt_Integration) # interval / timestep size?
  period= int(Iter_cycle/100)

  # Oscialltory
  if shear_style == 'oscillatory':
    #value, t_start, t_ramp
    shear_style_vinf = hoomd.variant.Cosinusoid(vinf, sim.timestep, t_ramp)
    shear_style_erate = hoomd.variant.Cosinusoid(shear_rate, sim.timestep, dt_Integration)

  # Sinusoid
  if shear_style == 'sinusoid': 
    #value, t_start, omega
    shear_style_vinf = hoomd.variant.Cosinusoid(vinf, sim.timestep, omega*dt_Integration)
    shear_style_erate = hoomd.variant.Cosinusoid(shear_rate, sim.timestep, omega*dt_Integration)

  # Cosinusoid
  if shear_style == 'cosinusoid' :
    #value, t_start, omega
    shear_style_vinf = hoomd.variant.Cosinusoid(vinf, sim.timestep, omega*dt_Integration)
    shear_style_erate = hoomd.variant.Cosinusoid(shear_rate, sim.timestep, omega*dt_Integration)


# set the random seed for reproducibility
seed_value = 42


######### SIMULATION
# Checks for existing shear flow files. If none exist, begins shearing 
# from the gelation state

if rank == 0 or comm.size == 1:
  if os.path.exists('Shear-'+shear_style+'-BD.gsd'):
    print('Shear flow file already exists. No new files created.')
    exit()
  else:
    print('Shearing the Brownian Dynamics simulation with '+shear_style+' shear')

## Create a CPU simulation
device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=seed_value)

# reset timestep to zero and start simulation from the GSD file
sim.timestep=0
sim.create_state_from_gsd(filename='Gelation-BD.gsd')

# assign particle types to groups
# (in case we want to integrate over subpopulations only,
# but would require other mods to source code)
groupA = hoomd.filter.Type(['B'])
all_ = hoomd.filter.Type(['B'])

# create neighboring list
nl = hoomd.md.nlist.Tree(buffer=0.05);

# define Morse force (attraction) interactions
# TO BE ADDED: bond_calc=bond_calc
morse = hoomd.md.pair.Morse(nlist=nl, default_r_cut=1.0 * r_c)

# colloid-colloid: hard particles (no deformation/overlap)
morse.params[('B', 'B')] = dict(D0=D0, alpha=kappa, r0=(R_C1+R_C1), f_contact=f_contact)
morse.r_cut[('A', 'A')] = (R_C1+R_C1) + (7.0/kappa) # used to assemble nl
#morse.r_cut[('B','B')] = r_c+(R_C1+R_C1) # used to assemble nl

# choose integration method for the end of each timestep
# BROWNIAN (overdamped) or LANGEVIN (underdamped)
brownian=hoomd.md.methods.Brownian(filter=all_, kT=kT, default_gamma=gamma)
integrator=hoomd.md.Integrator(dt=dt_Integration, vinf=shear_style_vinf, forces=[morse], methods=[brownian])
#langevin = hoomd.md.methods.Langevin(filter=all_, kT=KT, default_gamma=gamma)
#integrator = hoomd.md.Integrator(dt=dt_Integration, vinf=shear_style, forces=[morse], methods=[langevin])
sim.operations.integrator = integrator

# set the simulation to log certain values
logger = hoomd.logging.Logger()
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
sim.operations.computes.append(thermodynamic_properties)
logger.add(thermodynamic_properties,quantities=['kinetic_temperature','pressure_tensor','virial_ind_tensor','potential_energy'])
logger.add(sim,quantities=['tps'])

# set output file
gsd_writer = hoomd.write.GSD(trigger=t_ramp, filename='Shear-'+shear_style+'-BD.gsd', filter=all_, mode='wb', dynamic=['property','momentum','attribute'])
gsd_writer.write_diameter = True
#gsd_writer.maximum_write_buffer_size = 1e8 # max 100 million bytes
sim.operations.writers.append(gsd_writer)
gsd_writer.logger = logger

# set the box resize style
box_resize=hoomd.update.BoxShear(trigger=1, erate=shear_style_erate, deltaT=dt_Integration, flip=False)
sim.operations += box_resize

# run the simulation!
#sim.run(Iter_cycle*5)
n_cycles = 2
for i in range(n_cycles):
  sim.run(Iter_cycle)
  gsd_writer.flush()

if rank == 0 or comm.size == 1:	
  print('New '+shear_style+' shear Brownian Dynamics state (Shear-'+shear_style+'-BD.gsd) created.')

  # record sim parameters
  print("\nR_C1:", R_C1)
  #print("R_C2:", R_C2)
  #print("L_X:", L_X)
  print("phi:", phi)
  #print("percent_C1:", percent_C1)
  #print("percent_C2:", percent_C2)
  #print("poly:", poly)
  #print("poly_C1:", poly_C1)
  #print("poly_C2:", poly_C2)
  print("seed_value:", seed_value)
  print("D0/kT:", round(D0/KT))
  print("kappa:", kappa, "\n")
  print("L_Y:", L_Y)
  print("shear_style:", shear_style)
  if shear_style == "constant":
    print("vinf:", vinf)
    print("shear_rate:", shear_rate)
  else: 
    print("amplitude:", A)
    print("frequency:", omega)
    print("t_ramp:", t_ramp)
    print("vinf:", vinf)
    print("shear_rate:", shear_rate)


