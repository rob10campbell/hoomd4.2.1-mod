## Analyze the results of a BD shear simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 colloid type (typeid=0)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN CHECKS ON A SIMULATION section
"""
## This code performs the following shear specific analyses: 
##   - calculates:
##        * colloid velocity profile
##        * non-affinity (amount of affine vs. non-affine motion)
##        * the fabric tensor
##
## For other analyses, use the standard analysis module and analysis scripts
"""
## (Rob Campbell)


import numpy as np
import gsd.hoomd
import math
import module

##########################
""" INPUT PARAMETERS """
##########################
# path to sim data
filepath = '../Shear-BD.gsd'

# filepath to folder where data files will be created
data_outpath = 'data'

## For the noaffinity calculation
affine_ref_filepath = "../../3-gelation/Gelation-BD.gsd"
## calculate BDtimes for your sim
if shear_style = 'constant'
  dt_Integration = 1e-4
  theta = 1.0
  shear_rate = 0.1
  t_ramp = theta/shear_rate
else: 
  dt_Integration = 1e-3
  Amp = 0.1  
  omega = 0.01 
  shear_rate = Amp * omega
  t_ramp = Amp/shear_rate 
Iter_cycle = int(t_ramp / dt_Integration)
period = int(Iter_cycle/frames_per_strain) 
BDtime = period*dt_Integration


###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


#######
#Particle velocity profile
"""
# bin the velocities to create a velocity profile within your flow geometry
# NOTE: this calculates the COLLOID velocity profile
"""

# set the parameters for binning the system
nlayers_VP = 100 # n_layers in y direction (should be at least particle size)

def colloid_vel_profile_py(filename):
  traj = gsd.hoomd.open(filename, 'rb')
  nparticles=traj[-1].particles.N
  Lbox=traj[-1].configuration.box[0:3]
  vel=traj[-1].particles.velocity
  pos=traj[-1].particles.position
  typeid=traj[-1].particles.typeid
  module.vel_profile(nparticles,nlayers_VP,pos,vel,typeid,Lbox)
#######


#######
#Non-affinity
"""
# calculate if the motion/transformation is affine or non-affine
# by comparing the gelation with the last frame
"""
def non_affinity_py(filename, BDtime):
  # get the gel trajectory
  traj = gsd.hoomd.open(affine_ref_filepath, 'r')
  Lbox=traj[-1].configuration.box[:3]
  inv_Lbox = 1.0/Lbox
  index=[]
  index.extend(np.where(traj[-1].particles.typeid == [0])[0])
  init_pos=traj[-1].particles.position[index]
  pos_prev_unwrap=init_pos.copy()

  # get the shear trajectory
  traj = gsd.hoomd.open(filename, 'rb')
  nframe=len(traj)

  filename='non_affinity.txt'
  f=open(filename,'w')
  f.write("BDtime nonaff\n")

  # find how each particle moved
  for i in range(100):
    if(i==0):
      dr = traj[i].particles.position[index,:]-pos_prev_unwrap
    else:
      dr = traj[i].particles.position[index,:]-traj[i-1].particles.position[index,:]
    dr[:,2] -= Lbox[2]*np.rint(dr[:,2]*inv_Lbox[2])
    img = Lbox[1] * np.rint(dr[:,1]*inv_Lbox[1])
    dr[:,1] -= img
    dr[:,0] -= img*traj[i].configuration.box[3]
    dr[:,0] -= Lbox[0]*np.rint(dr[:,0]*inv_Lbox[0])
    pos_prev_unwrap += dr
    pos = pos_prev_unwrap - init_pos
    pos[:,0] -= i*BDtime*init_pos[:,1]

    # find if the motion is nonaffine
    nonaff = np.mean(np.linalg.norm(pos,axis=1))
    f.write('%f,%f\n'%(i*BDtime,nonaff))
  f.close()
#######


#######
# fabric tensor
"""
# calculate the fabric tenso
# see this paper for more info:
#     Jamali, McKinley, Armstrong (2017) PRL
#     "Microstructural Rearrangements and their Rheological Implications in a Model Thixotropic Elastoviscoplastic Fluid"
#     DOI: https://doi.org/10.1103/PhysRevLett.118.048003
"""
def fabric_tensor_py(filename):
  cut_off=0.3
  traj = gsd.hoomd.open(filename, 'r')
  nframe=len(traj)
  Lbox=traj[-1].configuration.box[:3]
  inv_Lbox = 1.0/Lbox
  index=[]
  index.extend(np.where(traj[-1].particles.typeid == [0])[0])
  nparticles=int(len(index))
  radius=0.5*traj[-1].particles.diameter[index]
  filename='fabric_tensor.txt'
  f=open(filename,'w')
  f.write("frame Rxy_avg\n")
  for i in range(0,nframe):
    pos=traj[i].particles.position[index]
    m_xy=traj[i].configuration.box[3]
    Rxy_avg=module.fabric_tensor(nparticles,pos,radius,Lbox,cut_off,m_xy)
    f.write("%f %f \n"%(i,Rxy_avg))
  f.close()
#######

if __name__ == '__main__':
  colloid_vel_profile_py(filename)
  non_affinity_py(filename,BDtime)
  fabric_tensor_py(filename)

