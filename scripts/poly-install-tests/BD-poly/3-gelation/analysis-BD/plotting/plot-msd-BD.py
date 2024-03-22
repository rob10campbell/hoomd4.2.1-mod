# plot the mean squared displacement (MSD) of colloids and solvent in
# a BD sim (from GSD file) versus time
#    (time can be frames, timesteps, or BD time)
#    options to include standard deviation or expected solvent MSD (y=x)
# NOTE: requires "msd.txt" created with the
#	sim-analysis-BD.py and module_analysis_BD.f90 scripts

import numpy as np
import math
import matplotlib.pyplot as plt

##############
""" INPUTS """
##############

# sim parameters
phi = 20
D0 = 12
R_C = 1
period = 10000
dt_Integration = 0.001
kT = 0.1

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: column names = frame msd_colloid msd_coll_smplstddev 
f = np.genfromtxt(data_directory+'/msd.txt', skip_header=1)

# to plot vs. number of timesteps: frame*period
timesteps = f[:,0]*period

# to plot vs. BD time: frame*t1
t1 = period * dt_Integration
BDtimes = f[:,0]*t1


def plot_msd_py(timechoice,flag,stddev):

  plt.plot(timechoice, f[:,1], 'grey', label="colloids")
  if stddev == True:
    plt.fill_between(timechoice, f[:,1]-f[:,2], f[:,1]+f[:,2],
      alpha=0.5, edgecolor='grey', facecolor='grey')

  # add particle size to the legend
  plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

  plt.title('Mean Squared Displacement', fontsize=16)
  if flag=='step':
    plt.xlabel('Number of timesteps', fontsize=16) 
  if flag=='BD':
    plt.xlabel('Brownian Dynamics (BD) times', fontsize=16)
  plt.ylabel('$\langle | x(t) - x_0 |^2 \\rangle$ [BD units]', fontsize=16)

  plt.loglog()

  # figsize includes title, axes, and plot; move lgd before this to include lgd in figsize 
  plt.rcParams['figure.figsize'] = [6, 6]
  plt.tight_layout()

  lgd = plt.legend(prop={"size":12}, loc='center left', bbox_to_anchor=(1,0.5), title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)

  plt.savefig('MSD_phi'+str(int(phi))+'_'+str(D0)+'kT.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600, transparent=False)
  #plt.show() #NOTE: plt.show() will cut off the legend, but the figure will save correctly

  print('MSD plot created')

if __name__ == '__main__':
   #plot_msd_py(timesteps,flag='step',stddev=False)
   plot_msd_py(BDtimes,flag='BD',stddev=False)
