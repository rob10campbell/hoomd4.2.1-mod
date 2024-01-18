# plot the average coordination number (<Z>) of colloids 
# in BD sim (from GSD file)
#    number of nearest neighbors versus time
#    (time can be frames, timesteps, or BD time)
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis-BD.py and module_analysis_BD.f90 scripts

import numpy as np
import pandas as pd
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

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = simframe, Z_any
df = pd.read_csv(data_directory+'/Zavg.csv')

frames = df['simframe'].to_numpy()

# to plot vs. number of timesteps: frame*period
timesteps = frames*period

# to plot vs. BD time: frame*t1
t1 = period * dt_Integration
BDtimes = frames*t1


def plot_Zavg_py(timechoice,flag):

   plt.plot(timechoice,df['Z_any'], 'grey', label="all colloids") # all colloids

   # add particle size to the legend
   plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

   #plt.ylim([0, 8])

   plt.title('Average Coordination Number', fontsize=16)
   if flag=='step':
      plt.xlabel('Number of timesteps', fontsize=16)          
   if flag=='BD':
      plt.xlabel('Brownian Dynamics (BD) times', fontsize=16)
   plt.ylabel('$\langle Z \\rangle$', fontsize=16)

   plt.legend(prop={"size":12}, loc=7, bbox_to_anchor=(1,0.4), title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
   plt.rcParams['figure.figsize'] = [6, 6]
   plt.tight_layout()

   plt.savefig('Zavg_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
   #plt.show()

   print('Zavg plot created')

if __name__ == '__main__':
   #plot_Zavg_py(timesteps,'step')
   plot_Zavg_py(BDtimes,'BD')
