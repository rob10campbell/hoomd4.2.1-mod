# plot the kinetic temperature (kT) of a simulation
# (from GSD file) versus time
#    (time can be frames, timesteps, or DPD time)
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis-DPD.py script

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
kT = 0.1 

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

# to plot vs. DPD time: frame*t1
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1


def plot_kT_py(timechoice,flag):

   # add phi and D0 label to legend
   plt.plot([], [], ' ', label='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT')

   plt.plot(timechoice,f[0:,8])

   plt.axhline(y=kT, color="r")
   plt.title('System Temperature', fontsize=16)
   if flag=='step':
      plt.xlabel('Number of timesteps', fontsize=16)
   if flag=='DPD':
      plt.xlabel('DPDtimes', fontsize=16)
   if flag=='diff':
      plt.xlabel('$t/\\tau_{C1}$ (colloid1 diffusion times)', fontsize=16)
   plt.ylabel('kT', fontsize=16)
 
   plt.legend(prop={"size":12})
   plt.rcParams['figure.figsize'] = [6, 6]
   plt.tight_layout()

   plt.savefig('kT_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
   #plt.show()

   print('Temperature (kT) plot created')

if __name__ == '__main__':
   #plot_kT_py(timesteps,'step')
   plot_kT_py(DPDtimes,'DPD')
