## Analyze the results of a DPD colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 solvent type (typeid=0) and 1 colloid type (typeid=1)
##
## Calculates:
## 	- [SHEARING] corrected (kinetic) pressure and temperature
##
## (Rob Campbell)


######### MODULE LIBRARY
import numpy as np
import gsd.hoomd
import math
import module

######### INPUT PARAMETERS

## corrected (kinetic) pressure and temperature
nlayerY = 50 # set umber of layers for binning/averaging in y-direction

######### DEFINE SIM CHECKS
## corrected (kinetic) pressure and temperature
## 1. extract ONLY the effect of colloid motion on pressure and temperature
## 2. normalize pressure by simulation volume

#  1. - divide the simulation into layers/bins along the y-axis (-L_Y/2 to L_Y/2)
#     - in each bin:
#         - find the average SOLVENT velocity in the x-direction
#         - subtract this bin/average solvent velocity from the x-velocity of EVERY 
#             particle in the bin (vel_x is now relative to solvent velocity)  
#         - use this corrected velocity to correct pressure (pressure is now ONLY
#	      the effect of colloid particles)
#     - use the corrected pressure to calculate kinetic temperature (kT)
#	(normalized by dimension (3) and number of particles)
#  2. - normalize the pressure by the simulation box volume


def corrected_temperature_py(filename):
  # open the simulations gsd file as "traj"
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframe=len(traj)
  # create a file to store corrected data
  f=open("T_corrected.txt",'w')
  # label each column (DPD timestep, pressure components, kinetic temperature)
  f.write("simframe pxx pxy pxz kT\n")

  # for each frame
  for i in range(0, nframe):
    # get the simulation parameters (set in init.gsd)
    nparticle = traj[i].particles.N # number of all particles
    typeid = traj[i].particles.typeid # typeids (solvent, colloid, etc)
    Lbox = traj[i].configuration.box[0:3] # simulation box size
    mass = traj[i].particles.mass # mass of each particle
    # and the current simulation state	
    pos = traj[i].particles.position # xyz position of each particle	
    vel = traj[i].particles.velocity # velocity vector of each particle

    # set the initial pressure components and temperature to zero
    pxx=0.0; pxy=0.0; pxz=0.0; kT=0.0

    # use the module to calculate corrected values with the assigned variables
    pxx,pxy,pxz,kT=module.corrected_temperature(nlayerY,nparticle,typeid,Lbox,mass,pos,vel)
    f.write("%f %f %f %f %f\n"%(i,pxx,pxy,pxz,kT))

  f.close()


######### RUN CHECKS ON A SIMULATION

if __name__ == '__main__':
	corrected_temperature_py('../../Shearing0.01-DPD.gsd')

