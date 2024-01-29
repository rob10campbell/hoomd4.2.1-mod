## serial simulation to create init.gsd
## for a DPD simulation with 1 solvent type and 1 colloid type
## NOTE: for a FIXED NUMBER OF COLLOIDS and variable box size
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
import random # psuedo-random number generator
# Other
import os # miscellaneous operating system interfaces


######### CUSTOM CLASSES
#N/A


#########  SIMULATION INPUTS
# General parameters
phi = 0.2 # volume fraction
rho = 3 # number density (per unit volume)

# Colloid particle details
N_C1 = 1000; # number of 1st type of colloidal particles
R_C1 = 1  # 1st type colloid particle radius
V_C1 = (4./3.) * math.pi * R_C1 ** 3 # 1st type colloid particle volume (1 particle)
m_C1 = V_C1 * rho # 1st type colloid particle mass

# Solvent particle details
m_S = 1 # solvent particle mass 
R_S = 0.5 # solvent particle radius 

# Simulation box size (calculated from # colloids)
L_X = (N_C1 * V_C1 / phi)**(1./3.)
L_Y = L_X
L_Z = L_X 

# Volumes
V_total = L_X * L_Y * L_Z # total volume of simulation box (cube)
# NOTE: volume calculation ASSUMES 1 colloid type and 1 solvent type
V_Colloids = N_C1 * V_C1 # total volume of all colloid particles
V_Solvents = V_total - V_Colloids # total volume of solvents 

# Calculate number of solvent particles (from volumes & number density)
N_Solvents = math.floor(rho * V_Solvents);

# Total number of particles in the simulation
N_total = int(N_Solvents + N_C1)

# set random seed for repeatable random number generation
seed_value = 42
numpy.random.seed(seed_value)

######### SIMULATION
## Checks for existing files. If none are found, creates a new 
## random distribution of particles for use as initial simulation state. 

if os.path.exists('init-DPD.gsd'):
  print("Initialization file already exists. No new files created.")
else:
  print("New DPD initialization file is being created")
  ## Initialize a snapshot of the system
  snapshot = gsd.hoomd.Frame()
  snapshot.configuration.box = [L_X, L_Y, L_Z, 0, 0, 0] # create the sim box
  snapshot.particles.N=N_total # add all particles to the snapshot
  # set the particle types
  snapshot.particles.types = ['A','B']
  # assign particles to each type
  typeid = []
  typeid.extend([0]*N_Solvents)
  typeid.extend([1]*N_C1)
  snapshot.particles.typeid = typeid
  # set a mass for each particle type
  mass = []
  mass.extend([m_S]*N_Solvents)
  mass.extend([m_C1]*N_C1)
  snapshot.particles.mass = mass
  # set a diameter for each particle type
  diameter = []
  diameter.extend([2.0*R_S]*N_Solvents)
  diameter.extend([2.0*R_C1]*N_C1)
  snapshot.particles.diameter = diameter
  # randomly distribute all the particles in 3D space
  pos_arr = numpy.zeros((N_total,3))
  pos_arr[:,0] = numpy.random.uniform(-0.5*L_X, 0.5*L_X, N_total)
  pos_arr[:,1] = numpy.random.uniform(-0.5*L_Y, 0.5*L_Y, N_total)
  pos_arr[:,2] = numpy.random.uniform(-0.5*L_Z, 0.5*L_Z, N_total)
  snapshot.particles.position = pos_arr
  # save the snapshot of the initialized system
  with gsd.hoomd.open(name='init.gsd', mode='w') as f:
    f.append(snapshot)

  print("New DPD initialization file (init-DPD.gsd) created.\n")
  print("Seed for Random Number Generator: "+str(seed_value)) 
  print("Simulation volume: L_X = " + str(L_X) + ", L_Y = " + str(L_Y) 
    + ", L_Z = " + str(L_Z))
  print("Volume fraction: " + str(phi) + "\n")
  print("Total number of particles: " + str(N_total))
  print("Total number of solvent particles: " + str(N_Solvents))
  print("Total number of colloid particles: " + str(N_C1))
