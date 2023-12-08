## mpi simulation to bring a DPD sim of an attractive 
## colloid system from equilibrium suspension to a 
## quasi-steady state gel
## NOTE: requires matching Equilibrium-DPD.gsd file
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

eta0 = 0.3 # background viscosity
gamma = 4.5 # DPD controlling parameter for viscous resistance (dissipative force)

N_time_steps = 1500000 # number of time steps
dt_Integration = 0.001 # dt! (DPD timestep)
period = 10000 # recording interval

# Colloid particle details
N_C1 = 1000 # number of 1st type of colloidal particles
R_C1 = 1  # 1st type colloid particle radius
V_C1 = (4./3.) * math.pi * R_C1 ** 3 # 1st type colloid particle volume (1 particle)
m_C1 = V_C1 * rho # 1st type colloid particle mass

# Solvent particle details
m_S = 1 # solvent particle mass 
R_S = 0.5 # solvent particle radius 

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  print('WARNING: r_c is less than range of attraction. Increase r_c')
r0 = 0.0 # minimum inter-particle distance
kappa = 30.0 # range of attraction (4 (long range)- 30 (short range)), distance in DPD units is approx 3/kappa
f_contact = 10000.0 * KT / r_c # set colloid-colloid hard-sphere interactions 
r_cut_sc = (r_c**3 + R_C1**3)**(1/3) # modified center-center cut-off radius for solvent-colloid interactions
bond_calc = False # do you want to track what bonds form and break? True=yes, False=no

# Simulation box size (NOTE: calculated from # colloids)
L_X = (N_C1 * V_C1 / phi)**(1./3.)
L_Y = L_X
L_Z = L_X 

# Volumes
V_total = L_X*L_Y*L_Z; # total volume of simulation box (cube)
# NOTE: volume calculation currently ASSUMES 1 colloid type and 1 solvent type
V_Colloids = N_C1 * V_C1; # total volume of all colloid particles
V_Solvents = V_total - V_Colloids; # total volume of solvents 

# Calculate number of solvent particles (from volumes & number density)
N_Solvents = math.floor(rho * V_Solvents);

# Total number of particles in the simulation
N_total = int(N_Solvents + N_C1)


######### SIMULATION
# Checks for existing gelation files. If none exist, continues from
# the equilibrium state towards a quasi-steady state gel (no shear).

if os.path.exists('Gelation-DPD.gsd'):
  print('Gelation file already exists. No new files created.')
else:
  print('DPD equilibrium state is being brought to quasi-steady state gelation')
  ## Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations

  # start the simulation from Equilibrium (and don't set timestep to 0)
  sim.create_state_from_gsd(filename='../2-equilibrium/Equilibrium-DPD.gsd')

  # assign particle types to groups
  # (in case we want to integrate over subpopulations only,
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['A'])
  groupB = hoomd.filter.Type(['B'])
  all_ = hoomd.filter.Type(['A','B'])

  # don't want to thermalize the system after equilibrium because
  # we are tracking velocity during gelation
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05);

  # define DPD Morse force (attraction) interactions
  morse = hoomd.md.pair.DPDMorse(nlist=nl, kT=KT, default_r_cut=1.0 * r_c, bond_calc=bond_calc)

  # solvent-solvent: soft particles (allow deformation/overlap)
  morse.params[('A','A')] = dict(A0=25.0 * KT / r_c, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=0.0, rcut=r_c) # force calc
  morse.r_cut[('A','A')] = r_c # used to assemble nl

  # solvent-colloid: soft particles (allow deformation/overlap)
  morse.params[('A','B')] = dict(A0=25.0 * KT / r_cut_sc, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, a1=0.0, a2=R_C1, 
    rcut=r_cut_sc - (0 + R_C1)) # force calc
  morse.r_cut[('A','B')] = r_cut_sc # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','B')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C1, a2=R_C1, rcut=r_c) # force calc
  morse.r_cut[('B','B')] = (r_c + 2.0 * R_C1) # used to assemble nl

  # choose integration method for the end of each timestep
  nve = hoomd.md.methods.ConstantVolume(filter=all_, thermostat=False)
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
  gsd_writer = hoomd.write.GSD(trigger=period, filename="Gelation-DPD.gsd", 
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
	
  print('New DPD gelation state (Gelation-DPD.gsd) created.')
