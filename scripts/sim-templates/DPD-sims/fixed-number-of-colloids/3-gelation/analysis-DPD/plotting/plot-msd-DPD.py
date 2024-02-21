# plot the mean squared displacement (MSD) of colloids and solvent in
# a DPD sim (from GSD file) versus time
#    (time can be frames, timesteps, or DPD time)
#    options to include standard deviation or expected solvent MSD (y=x)
# NOTE: requires "msd.txt" created with the
#	sim-analysis-DPD.py and module_analysis_DPD.f90 scripts

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

# load data: column names =
# frame msd_colloid msd_coll_smplstddev msd_solvent msd_solv_smplstddev
f = np.genfromtxt(data_directory+'/msd.txt', skip_header=1)

# to plot vs. number of timesteps: frame*period
timesteps = f[:,0]*period

# to plot vs. DPD time: frame*t1
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1


def plot_msd_py(timechoice,flag,solvent,stddev,expected_solvent):

  plt.plot(timechoice, f[:,1], 'grey', label="colloids")
  if stddev == True:
    plt.fill_between(timechoice, f[:,1]-f[:,2], f[:,1]+f[:,2],
      alpha=0.5, edgecolor='grey', facecolor='grey')

  if solvent == True:
    plt.plot(timechoice, f[:,3], 'blue', label="solvent")
    if stddev == True:
      plt.fill_between(timechoice, f[:,3]-f[:,4], f[:3]+f[:,4],
        alpha=0.5, edgecolor='blue', facecolor='blue')

  # compare with expected solvent
  if expected_solvent == True:
    # plot y=x (the expected slope of solvent MSD)
    x = np.linspace(10, 1990)
    plt.plot(x, x, linestyle='--', color='lightgrey')

  # add particle size to the legend
  plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

  plt.title('Mean Squared Displacement', fontsize=16)
  if flag=='step':
    plt.xlabel('Number of timesteps', fontsize=16) 
  if flag=='DPD':
    plt.xlabel('DPD times', fontsize=16)
  plt.ylabel('$\langle | x(t) - x_0 |^2 \\rangle$ [DPD units]', fontsize=16)

  plt.loglog()

  plt.legend(prop={"size":12}, loc=7, bbox_to_anchor=(1,0.4), title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
  plt.rcParams['figure.figsize'] = [6, 6]
  plt.tight_layout()

  plt.savefig('MSD_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
  #plt.show()

  print('MSD plot created')

if __name__ == '__main__':
   #plot_msd_py(timesteps,flag='step',solvent=False,stddev=False,expected_solvent=False)
   plot_msd_py(DPDtimes,flag='DPD',solvent=False,stddev=False,expected_solvent=False)
