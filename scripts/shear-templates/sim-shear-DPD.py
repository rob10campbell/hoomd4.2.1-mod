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

# use higher viscosity to see shear behavior better
eta0 = 1.1 # background viscosity
gamma = 45 # DPD controlling parameter for viscous resistance (dissipative force)

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if(snap.communicator.rank == 0): # if this is the start of the simulation
  if r_c < (3/kappa):
    print('WARNING: r_c is less than range of attraction. Increase r_c')
r0 = 0.0 # minimum inter-particle distance
f_contact = 10000.0 # magnitude of contact force (usually 100 or 1000)
bond_calc = False # do you want to track what bonds form and break? True=yes, False=no

# modified center-center cut-off radius for solvent-colloid interactions
r_sc1_cut = (r_c**3 + R_C1**3) ** (1/3)
r_sc2_cut = (r_c**3 + R_C2**3) ** (1/3)

# shear flow details
shear_style = 'constant'
#shear_style = 'cosinusoid' 

n_strains = 1 # number of strains or oscillation cycles
             # ex: 10 for constant strain, 5 for oscillatory, etc.
frames_per_strain = 100 # how many frames you want in each cycle
n_frames = frames_per_strain*n_strains # total simulation frames

init_velocity = False # initialize with or without a linear velocity profile

# box size
#L_X = 30 # box size in flow direction
L_Y = 30 # box size in gradient direction
#L_Z = 30 # box size in vorticity direction

# Constant
if shear_style == 'constant':
  # NOTE: reduce timestep 0.001->0.000001 as SR increases 0.001->1.0
  dt_Integration = 1e-4 # dt! (BD timestep, may need to be smaller than DPD)

  shear_rate = 0.1
  vinf = shear_rate * L_Y # the maximum flow velocity
  theta = 1.0 # desired xy tile factor (0.0 to 1.0)
  t_ramp = theta/shear_rate # length of each shear pass in BD-time

else:
  # NOTE: reduce timestep 0.001->0.000001 as SR increases 0.001->1.0
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

if(snap.communicator.rank == 0): # if this is the start of the simulation
  if os.path.exists('Shear-DPD.gsd'):
    print('Shear flow file already exists. No new files created.')
    exit()
  else:
    print('Shearing the DPD simulation with '+shear_style+' shear')

## Create a CPU simulation
device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=seed_value)

# reset timestep to zero and start simulation from the GSD file
sim.timestep=0
sim.create_state_from_gsd(filename='../3-gelation/Gelation-DPD.gsd')

# assign particle types to groups
# (in case we want to integrate over subpopulations only,
# but would require other mods to source code)
groupA = hoomd.filter.Type(['A'])
groupB = hoomd.filter.Type(['B'])
groupC = hoomd.filter.Type(['C'])
all_colloids = hoomd.filter.Type(['B','C'])
all_ = hoomd.filter.Type(['A','B','C'])

# DON'T thermalize the system (shearing is not at thermal equilibrium!)
#sim.state.thermalize_particle_momenta(filter=groupA, kT=kT)

# [optional] apply initial velocity field
if init_velocity == True:
  # apply an initial velocity field to all solvent particles
  snap = sim.state.get_snapshot()
  if(snap.communicator.rank == 0): # if this is the start of the simulation
    for i in range(snap.particles.N): # for all particles
      #if(snap.particles.typeid[i]==0): # if that particle is a solvent particle
      v_x = snap.particles.position[i,1]*SR # apply a linear velocity profile
      snap.particles.velocity[i][0] += v_x
  # apply the linear velocity profile to the initial simulation state
  sim.state.set_snapshot(snap)

# set shear style Variant class
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
morse = hoomd.md.pair.DPDMorse(nlist=nl, kT=kT, default_r_cut=1.0 * r_c, bond_calc=bond_calc)

# solvent-solvent: soft particles (allow deformation/overlap)
morse.params[('A','A')] = dict(A0=25.0 * kT / r_c, gamma=gamma,
  D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0,
  a1=0.0, a2=0.0, rcut=r_c) # force calc
#morse.r_cut[('A', 'A')] = r_c + (7.0/kappa) # used to assemble nl
morse.r_cut[('A','A')] = r_c # used to assemble nl

# solvent-colloid: soft particles (allow deformation/overlap)
morse.params[('A','B')] = dict(A0=25.0 * kT / r_sc1_cut, gamma=gamma,
  D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, a1=0.0, a2=R_C1,
  rcut=r_sc1_cut - (0 + R_C1)) # force calc
#morse.r_cut[('A', 'B')] = r_sc1_cut + (7.0/kappa) # used to assemble nl
morse.r_cut[('A','B')] = r_sc1_cut # used to assemble nl

# colloid-colloid: hard particles (no deformation/overlap)
morse.params[('B','B')] = dict(A0=0.0, gamma=gamma,
  D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact,
  a1=R_C1, a2=R_C1, rcut=r_c) # force calc
#morse.r_cut[('B', 'B')] = (r_c + R_C1+R_C1) + (7.0/kappa) # used to assemble nl
morse.r_cut[('B','B')] = (r_c + 2.0 * R_C1) # used to assemble nl

# colloid-colloid: hard particles (no deformation/overlap)
morse.params[('C','C')] = dict(A0=0.0, gamma=gamma,
  D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact,
  a1=R_C2, a2=R_C2, rcut=r_c) # force calc
#morse.r_cut[('C', 'C')] = (r_c + R_C2+R_C2) + (7.0/kappa) # used to assemble nl
morse.r_cut[('C','C')] = (r_c + 2.0 * R_C2) # used to assemble nl

# choose integration method for the end of each timestep
nve = hoomd.md.methods.ConstantVolume(filter=all_, thermostat=False)
integrator=hoomd.md.Integrator(dt=dt_Integration, vinf=shear_style_vinf, forces=[morse], methods=[nve])
sim.operations.integrator = integrator

# set the simulation to log certain values
logger = hoomd.logging.Logger()
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
sim.operations.computes.append(thermodynamic_properties)
logger.add(thermodynamic_properties,quantities=['kinetic_temperature','pressure_tensor','virial_ind_tensor','potential_energy'])
logger.add(sim,quantities=['tps'])

# set output file
gsd_writer = hoomd.write.GSD(trigger=period, filename='Shear-DPD.gsd', filter=all_, mode='wb', dynamic=['property','momentum','attribute'])
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
  print('\nNew '+shear_style+' shear DPD state (Shear-DPD.gsd) created.')

  if shear_style == 'constant':
   print(shear_style+' shear')
   print('initial velocity? '+str(init_velocity))
   print('shear rate, SR: '+str(round(shear_rate,2))
   print('t_ramp:', str(t_ramp))
   print('dt_Integration:', str(dt_Integration))
   print('frames per strain:', str(frames_per_strain))
   print('Iter_cycle:', str(Iter_cycle))
   print('period:', str(period))
   print('n_strains:',str(n_strains))

  else:
   print(shear_style+' shear')
   print('initial velocity? '+str(init_velocity))
   print('Amplitude: '+str(Amp,3))
   print('omega: '+str(omega,3))
   print('shear rate, SR: '+str(round(shear_rate,2))
   print('t_ramp:', str(t_ramp))
   print('dt_Integration:', str(dt_Integration))
   print('frames per strain:', str(frames_per_strain))
   print('Iter_cycle:', str(Iter_cycle))
   print('period:', str(period))
   print('n_strains:',str(n_strains))

