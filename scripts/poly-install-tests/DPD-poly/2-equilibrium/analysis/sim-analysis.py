## Analyze the results of a DPD colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 solvent type (typeid=0) and 1 colloid type (typeid=1)
##
## Extracts temperature and pressure (-stress) data from a GSD file, 
## then calculates:
## 	- average coordination number (Z) of a colloid per simulation frame
##	- mean squared displacement (MSD) of all colloids in each frame
##	- pair correlation function (PCF) g(r1, r2): the probability density per unit 
##		volume of a particle at r1 finding another particle at r2
##	- radial distribution fuction (RDF) g(r): the density of (other) particles per 
##		unit volume relative to the position of a chosen particle (r1)
##
## (Rob Campbell)


######### MODULE LIBRARY
import numpy as np
import gsd.hoomd
import math
import module

######### INPUT PARAMETERS
filepath = '../Equilibrium.gsd'

## general simulation parameters
period = 10000 # from simulation
dt_Integration = 0.001 # from simulation
t1 = period * dt_Integration # timestep conversion factor

## coordination number parameters
cut_off = 0.1 # cut-off radius for Z calculation

## msd parameters
Lbox_shortest = 30 # from simulation (shortest side)
kT = 0.1 # from simulation
R_C1 = 1 # from simulation
eta0 = 1.1 # from simulation
d = 3 # dimension of the system (2D = 2, 3D = 3)
D = kT/(6*math.pi*eta0*R_C1) # diffusion coefficient = r^2 / 2d*tau
tau_to_half = (Lbox_shortest/2)**2 / (2*d*D) # diffusion time to half-box (L/2)

## PCF parameters
# TBD

## RDF parameters
Lbox = 30 # from simulation (assumes cube)
V_total = Lbox**3 # total volume (assumes cube)
rho = 3 # from simulation
rdf_rmax = 6.0
rdf_bin_width = 0.1


######### DEFINE SIM CHECKS
## extract thermodynamic properties (temperature and pressure (AKA -stress) components)
def extract_properties_py(filename):
	# open the simulation GSD file as "traj" (trajectory)
	traj = gsd.hoomd.open(filename, 'r')
	# get the number of frames
	nframe = len(traj)
	
	# create a file to save the thermodynamic properties
	f=open('gsd-properties.txt', 'w')
	f.write('DPD-time Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps\n')

	# for each frame
	for i in range(0, nframe):
		# extract the total "force virial contribution" of the pressure tensor
		Pr=float(traj[i].log['md/compute/ThermodynamicQuantities/pressure_tensor'][1])
		# extract the decomposed virial compoenents of the pressure tensor
		Vr=np.asarray(traj[i].log['md/compute/ThermodynamicQuantities/virial_ind_tensor'],dtype=float)
		# extrac the potential energy and scale to kilo-units (1/1000)
		Pe=float(traj[i].log['md/compute/ThermodynamicQuantities/potential_energy'][0])*0.001
		# extract the kinetic temperature (kT)
		KT=float(traj[i].log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0])
		# extract the transactios per secont (tps) for assessing speed/efficiency of the simulation
		tps=int(traj[i].log['Simulation/tps'][0])

		# write these values to the output file:
		# raw values
		f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format((i+1)*t1, Pr, 
			Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))
		
		# rounded values
		#f.write('%f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %d %0.2f %d\n'%((i+1)*t1, Pr, 
		#	Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))


## average coordination number (Z) per frame
# calclates the average coordination number (Z) of a colloid particle in each frame
# of a simulation GSD file
def coordination_number_py(filename):
	# open the simulation GSD file as "traj" (trajectory)
	traj = gsd.hoomd.open(filename, 'r')
	# get the number of frames
	nframe = len(traj)
	# get the simulation box size [L_X, L_Y, L_Z] from the last frame 
	Lbox = traj[-1].configuration.box[:3]
	# get all the colloid particles in the last frame
	colloids = np.where(traj[-1].particles.typeid == [1])[0]
	# use this to count the total number of colloids
	ncolloid = len(colloids)
	# use last frame to calculate radius of each colloid
	R_C1 = 0.5*traj[-1].particles.diameter[colloids]

	# create a file to hold the coordination number data
	f=open('Zavg.txt','w')
	f.write('DPD-time Zavg\n')

	# for each frame
	for i in range(0,nframe):
		# find all the colloid particles
		colloids = np.where(traj[i].particles.typeid == [1])[0]
		# get their positions
		pos = traj[i].particles.position[colloids]
		# get the xy tilt factor (square=0.0, sheared-right=0.45)
		m_xy = traj[i].configuration.box[3]
		# run the module to calculate the Zavg for the frame
		Zavg = module.coordination_number(Lbox,ncolloid,R_C1,pos,m_xy,cut_off)

		# write the data to the output file
		f.write("{0} {1}\n".format(round((i+1)*t1),Zavg))
	f.close()


## mean squared displacement (msd)
# calculates the mean squared displacement (MSD) for all colloid particles in each frame
# of a simulation GSD file, as well as the (corrected) sample standard deviation of the MSD
def msd_py(filename):
	if period >= tau_to_half:
		print('Error: the GSD file\'s recording timestep is too large for MSD calculations.\n\n' + 
			'To create data that you can use to calculate the MSD, you should set the'
			' period/trigger to a small enough value that a particle should not move'
			' 1/2 the box length in between frames.\n\nRerun the simulation with' +
			' a smaller trigger/period before calculating MSD.')
	else:
		# open the simulation GSD file as "traj" (trajectory)
		traj = gsd.hoomd.open(filename, 'r')
		# get the simulation box size [L_X, L_Y, L_Z] from the last frame 
		Lbox = traj[-1].configuration.box[:3]
		inv_Lbox = 1.0/Lbox
		# get the number of frames
		nframe = len(traj)	
		# get all the colloid particles in the last frame
		colloids = np.where(traj[-1].particles.typeid == [1])[0]
		# use this to count the total number of colloids
		ncolloid = len(colloids)

		# create an empty array for xyz positon of all colloids in all frames 	
		allpos = np.zeros((nframe,ncolloid,3))
		# get the initial colloid positions from the first frame
		allpos[0,:,:] = traj[0].particles.position[colloids] 
		# correct the change in position for colloids crossing a box boundary 
		for i in range(1,nframe):
			# calculate the change in position since the previous frame
			delpos = traj[i].particles.position[colloids] - traj[i-1].particles.position[colloids]
			# if it is more than one box length (i.e. the particle crossed a 
			# boundary) correct value the position to be inside the sim box
			delpos -= Lbox*np.rint(delpos*inv_Lbox)
			# update and record the position of the colloids in this frame
			realpos = allpos[i-1,:,:] + delpos
			allpos[i,:,:] = realpos

		# run the module to calculate MSD
		# (output file "msd.txt" is created in the Fortran module)	
		module.msd_calculation(nframe,ncolloid,allpos,t1) 


## pair correlation function (PCF) g(r1, r2)
# the probability density per unit volume of a particle at r1 finding another particle at r2
# works for all systems, depends on both of the particles' positions (i.e. the magnitude AND 
# direction of the interparticle separation vectory)
#
#   PCF g(r1, r2) returns a probability density (per unit volume) normalized by 
#   the number density at both r1 and r2 (r1 and r2 are both position VECTORS)




## radial distribution function (RDF) g(r)
## the relative density of (other) particles per unit volume relative (normalized by the
## density of an isotropic system, number-density*volume) at a relative position away from a 
## chosen particle (r1)
# works for homogeneous, translationally invariant, isotropic systems, depends on separation 
# distance (i.e. only the magnitude, NOT direction, of the interparticle separation vector)
# for homogeneous, translationally invariant systems:
#   g(r) ≡ g(r1-r2) (r is a VECTOR relative to r1; r1 and r2 are position VECTORS)
#   number density ⟨n(r)⟩ = ⟨n⟩ (is independent of the vector r)
# for isotropic systems: g(r) is a function of the separation distance (a SCALAR: r =|r|)
#
#   RDF g(r) returns a density per unit volume at position r (r is a SCALAR relative to position r1)
#		normalized by the density (in that volume) of an isotropic system
def rdf_py(filename):	
	# open the simulation GSD file as "traj" (trajectory)
	traj = gsd.hoomd.open(filename, 'r')
	# get the simulation box size [L_X, L_Y, L_Z] from the last frame 
	Lbox = traj[-1].configuration.box[:3]
	# get the xy tilt factor (box deformation in x direction) from the last frame
	m_xy = traj[-1].configuration.box[3]
	# get the total number of particles from the last frame
	nparticle = traj[-1].particles.N
	# get all the colloid particles in the last frame
	colloids = np.where(traj[-1].particles.typeid == [1])[0]
	# use this to count the total number of colloids
	ncolloid = len(colloids)
	# number density of solvent (from the simulation)
	rho_S = rho
	# number density of colloid
	rho_C = ncolloid / V_total

	# get the typeid of all particles from the last frame
	typeid = traj[-1].particles.typeid
	# get the xyz position of all particles from the last frame
	pos=traj[-1].particles.position
	# set the maximum distance considered from the center of a particle
	# (AKA size of the search box)
	rmax = rdf_rmax
	# set the size of each bin/layer/slab
	bin_width = rdf_bin_width
	# use rmax and bin_width to calculate the number of layers
	nlayers = int(round(rmax/bin_width))

	# run the module to calculate g(r)
	# (output files "rdf-CS.txt" and "rdf-CC" are created in the Fortran module)	
	module.rdf_calc(Lbox,m_xy,nparticle,ncolloid,rho_S,rho_C,typeid,pos,bin_width,nlayers)
d

######### RUN CHECKS ON A SIMULATION

if __name__ == '__main__':	
	extract_properties_py(filepath)
	#coordination_number_py(filepath)
	#msd_py(filepath)
	#rdf_py(filepath)
