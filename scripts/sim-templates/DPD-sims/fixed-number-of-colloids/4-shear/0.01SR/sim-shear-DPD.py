## mpi simulation to apply shear flow to
## a DPD simulation of a quasi-steady state gel
## NOTE: requires matching Gelation-DPD.gsd file
##       for a FIXED NUMBER OF COLLOIDS and variable box size
## (Rob Campbell)

######### MODULE LIBRARY
# Use HOOMD-blue
import hoomd
import hoomd.md # molecular dynamics
import hoomd.custom
# Use GSD files
import gsd # provides Python API; MUST import sub packages explicitly (as below)
import gsd.pygsd # GSD reader in pure Python
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy
import math
import random # psuedo-random number generator
# Other
import os # miscellaneous operating system interfaces


######### CUSTOM CLASSES
# create a class to rotate the simulation box 180 degrees
# (after it has been deformed: /_/ -> \_\ )
## this prevents bonds from artifically breaking and reforming at the boundary
## and seems to provide correct shear properties
## HOWEVER, it switches to oscillatory flow.... might destroy strucure info?
#class flip_box(hoomd.custom.Action):
#	def act(self, timestep):
#		# get a snapshot of the current system /_/
#		snap = self._state.get_snapshot()
#		# if we're at the start of a shearing step
#		if(snap.communicator.rank == 0):
#			# select the simulation box
#			box_ = snap.configuration.box
#			# flip the sign of the xy tilt factor (tilt the box back to neutral)
#			snap.configuration.box = [box_[0], box_[1], box_[2], 
#							-box_[3], box_[4], box_[5]]
#			# get the original particle positions
#			x = snap.particles.position[:,0]
#			z = snap.particles.position[:,2]
#			# update positions to the new coordinate system
#			snap.particles.position[:,0] = -x
#			snap.particles.position[:,2] = -z
#		# update the simulation
#		self._state.set_snapshot(snap)


######### SIMULATION INPUTS 
# General parameters
phi = 0.2; # volume fraction
rho = 3.0 # number density (per unit volume)
KT = 0.1 # system temperature
D0 = 12.0 * KT # attraction strength (gels at >=4kT)

# use higher viscosity to see shear behavior better
eta0 = 1.1 # background viscosity
gamma = 45 # DPD controlling parameter for viscous resistance (dissipative force)

# N_time_steps //calculated later based on shearing parameters//
dt_Integration = 0.001; # dt! (DPD timestep)
period = 1 # recording interval

SR = 0.01 # shear parameter (shear rate = SR*L_X)
N_strains = 10 # number of strains

# set the remaining shearing parameters
theta = 1.0 # xy tilt factor
delta_T_shearing = round(theta/SR/dt_Integration) # timestep for shear
nframe_strain = delta_T_shearing/10 # 20 frames in theta strains

# Colloid particle details
N_C1 = 1000 # number of 1st type of colloidal particles
R_C1 = 1.0; # 1st type colloid particle radius
V_C1 = (4./3.) * math.pi * R_C1 ** 3; # 1st type colloid particle volume (1 particle)
m_C1 = V_C1 * rho; # 1st type colloid particle mass

# Solvent particle details
m_S = 1.0; # solvent particle mass
R_S = 0.5; # solvent particle radius

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter (r_cut = # * r_c)
r0 = 0.0 * R_C1 # minimum inter-particle distance 
kappa = 30. # range of attraction (4 (long range)- 30 (short range)), distance in DPD units is approx 3/kappa 
f_contact = 10000. * KT / r_c # set colloid-colloid hard-sphere interactions 
r_cut_sc = (r_c**3 + R_C1**3)**(1/3) # modified center-center cut-off radius for solvent-colloid interactions
bond_calc = False # do you want to track what bonds form and break? True=yes, False=no

# Simulation box size (NOTE: calculated from # colloids)
L_X = (N_C1 * V_C1 / phi)**(1./3.)
L_Y = L_X
L_Z = L_X

# Volumes
V_total = L_X*L_Y*L_Z # total volume of simulation box (cube)
# NOTE: volume calculation currently ASSUMES 1 colloid type and 1 solvent type
V_Colloids = N_C1 * V_C1 # total volume of all colloid particles
V_Solvents = V_total - V_Colloids # total volume of solvents (for 1 colloid type and 1 solvent type)

# Calculate number of solvent particles (from volumes & number density)
N_Solvents = math.floor(rho * V_Solvents);  

# Total number of particles in the simulation
N_total = int(N_Solvents + N_C1)


######### SIMULATION
# Checks for exisiting shearing files at the current shear rate (SR). If none exist,
# applies a velocity equal to the shear rate at the +Y and -Y boundaries and fully 
# deforms the system a set number of times (N_strains) to shear the gel.

if os.path.exists('Shearing' + str(SR) + '-DPD.gsd'):
  print('Shearing file already exists for shear rate ' + str(SR) + '. No new files created')
else:
  print('Shearing the DPD system at ' + str(SR))
  # Create a CPU simulation
  device = hoomd.device.CPU()
  sim = hoomd.Simulation(device=device, seed=50) # set seed to a fixed value for reproducible simulations
	
  # OPTIONAL? set the initial timestep to zero
  #sim.timestep=0
	
  # start shearing from the quasi-steady state gel created in the gelation step
  sim.create_state_from_gsd(filename='../../3-gelation/Gelation-DPD.gsd')

  # assign particle types to groups
  # (in caswe we want to integrate over subpopulations only,
  # but would require other mods to source code)
  groupA = hoomd.filter.Type(['A'])
  groupB = hoomd.filter.Type(['B'])
  all_ = hoomd.filter.Type(['A','B'])

  # DON'T thermalize the system (shearing is not at thermal equilibrium!)
  #sim.state.thermalize_particle_momenta(filter=all_, kT=KT)

  # apply an initial velocity field to all solvent particles
  snap = sim.state.get_snapshot()
  if(snap.communicator.rank == 0): # if this is the start of the simulation
    for i in range(snap.particles.N): # for all particles
      #if(snap.particles.typeid[i]==0): # if that particle is a solvent particle
      v_x = snap.particles.position[i,1]*SR # apply a linear velocity profile
      snap.particles.velocity[i][0] += v_x
  # apply the linear velocity profile to the initial simulation state
  sim.state.set_snapshot(snap)

  # create neighboring list
  nl = hoomd.md.nlist.Tree(buffer=0.05)

  # define DPD Morse force (attraction) interactions
  morse = hoomd.md.pair.DPDMorse(nlist=nl, kT=KT, default_r_cut=1.0 * r_c, bond_calc=bond_calc)

  # solvent-solvent: soft particles (allow deformation/overlap)
  morse.params[('A','A')] = dict(A0=25.0 * KT / r_c, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=0.0, rcut=r_c) # force calc
  morse.r_cut[('A','A')] = r_c # used to assemble nl

  # solvent-colloid: soft particles (allow deformation/overlap)
  morse.params[('A','B')] = dict(A0=25.0 * KT / r_cut_sc, gamma=gamma, 
    D0=0, alpha=kappa, r0=r0, eta=0.0, f_contact=0.0, 
    a1=0.0, a2=R_C1, rcut=r_cut_sc - (0 + R_C1)) # force calc
  morse.r_cut[('A','B')] = r_cut_sc # used to assemble nl

  # colloid-colloid: hard particles (no deformation/overlap)
  morse.params[('B','B')] = dict(A0=0.0, gamma=gamma, 
    D0=D0, alpha=kappa, r0=r0, eta=eta0, f_contact=f_contact, 
    a1=R_C1, a2=R_C1, rcut=r_c) # force calc
  morse.r_cut[('B','B')] = (r_c + 2.0 * R_C1) # used to assemble nl

  # choose integration method at the end of each timestep
  nve = hoomd.md.methods.ConstantVolume(filter=all_, thermostat=None)
  integrator = hoomd.md.Integrator(dt=dt_Integration, SR=SR*L_Y, forces=[morse], methods=[nve])
  sim.operations.integrator = integrator

  # set the simulation to log additional values in the gsd file
  logger = hoomd.logging.Logger()
  thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
  sim.operations.computes.append(thermodynamic_properties)
  logger.add(thermodynamic_properties, quantities=['kinetic_temperature',
    'pressure_tensor','virial_ind_tensor','potential_energy'])
  logger.add(sim, quantities=['tps'])

  # save outputs
  gsd_writer = hoomd.write.GSD(trigger=nframe_strain, filename='Shearing'+str(SR)+'-DPD.gsd', 
  filter=all_, mode='wb', dynamic=['property','momentum','attribute']) 
  gsd_writer.write_diameter = True
  sim.operations.writers.append(gsd_writer)
  gsd_writer.logger = logger

  ## set the box resize operation
  initial_box = hoomd.Box.from_box(sim.state.box)
  final_box = hoomd.Box.from_box(initial_box)
  final_box.xy = theta
  delta_T_shearing = round((final_box.xy-initial_box.xy)/SR/dt_Integration)
  box_resize=hoomd.update.BoxResize(trigger=period, box1=initial_box, box2=final_box,
    variant=hoomd.variant.Ramp(A=0, B=1, t_start=sim.timestep, t_ramp=delta_T_shearing),
    filter=all_, SR=SR*L_Y)

  ## set the box flip operation if using
  #flip_operation = hoomd.update.CustomUpdater(action=flip_box(), 
  #  trigger=hoomd.trigger.On(sim.timestep+delta_T_shearing))


  # append these operations to the system
  sim.operations += box_resize
  #sim.operations += flip_operation

  # shear the box one step (theta = 1.0)
  sim.run(delta_T_shearing+1)

  # write data from this step to the gsd file
  gsd_writer.flush()

  # clear the operations for the next round
  sim.operations -= box_resize
  #sim.operations -= flip_operation

  # shear the system
  for i in range(N_strains):
    # set the box resize operation
    initial_box = hoomd.Box.from_box(sim.state.box)
    initial_box.xy = -theta # reset the box to zero
    final_box = hoomd.Box.from_box(initial_box)
    final_box.xy = theta
    delta_T_shearing = round((final_box.xy-initial_box.xy)/SR/dt_Integration)
    box_resize = hoomd.update.BoxResize(trigger=period, box1=initial_box, box2=final_box,
      variant = hoomd.variant.Ramp(A=0, B=1, t_start=sim.timestep, 
      t_ramp=delta_T_shearing), filter=all_, SR=SR*L_Y)

    # set the flip operations if using
    #flip_operation = hoomd.update.CustomUpdater(action=flip_box(), 
    #  trigger=hoomd.trigger.On(sim.timestep+delta_T_shearing))

    # add the operations to the system
    sim.operations += box_resize	
    #sim.operations += flip_operation

    # run for one shear times (theta = 1.0, for a full 1 strain)
    sim.run(delta_T_shearing+1)

    # write data from this step to the gsd file
    gsd_writer.flush()

    # clear the operations for the next strain
    sim.operations -= box_resize
    #sim.operations -= flip_operation

    print('DPD shearing data saved to Shearing' + str(SR) + '-DPD.gsd')

