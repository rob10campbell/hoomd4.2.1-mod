# plot the average shear stress of all particles (from GSD file) versus time
#    (time can be frames, timesteps, or BD time)
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis-BD.py script

import numpy as np
import math
import matplotlib.pyplot as plt

##############
""" INPUTS """
##############

# sim parameters
phi = 20
D0 = 12
period = 10000
dt_Integration = 0.001

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: column names =
# simframe Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps
f = np.genfromtxt(data_directory+'/gsd-properties.txt', skip_header=1)

#to plot t vs. number of timesteps: frame*period
timesteps = f[:,0]*period

# to plot vs. BD time: frame*t1
t1 = period * dt_Integration
BDtimes = f[:,0]*t1

def plot_shearstress_py(timechoice,flag):

   # add phi and D0 label to legend
   plt.plot([], [], ' ', label='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT')

   plt.plot(timechoice,-f[0:,1])

   plt.axhline(y=0.0, color="r")
   plt.title('System Shear Stress', fontsize=16)
   if flag=='step':
      plt.xlabel('Number of timesteps', fontsize=16)
   if flag=='BD':
      plt.xlabel('Brownian Dynamics (BD) times', fontsize=16)
   plt.ylabel('$\sigma$', fontsize=16)

   plt.legend(prop={"size":12})
   plt.tight_layout()

   plt.savefig('shearstress_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
   #plt.show()

   print('Shear stress plot created')

if __name__ == '__main__':
   #plot_shearstress_py(timesteps,'step')
   plot_shearstress_py(BDtimes,'BD')
