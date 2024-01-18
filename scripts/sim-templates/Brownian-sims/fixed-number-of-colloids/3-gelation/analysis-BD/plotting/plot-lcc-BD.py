# plot the largest connected compoinent (LCC) of colloids 
# in BD sim (from networkX analysis)
#    percent of colloids in the largest connected component versus time
#    (time can be frames, timesteps, or BD time)
# NOTE: requires "networkx-allframes.csv" created with the
#	sim-networkCSV-BD.py script

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import os

##############
""" INPUTS """
##############

# sim parameters
phi = 20
D0 = 12
R_C = 1
period = 10000
dt_Integration = 0.001

# choose timeaxis for plot:  'frame', 'timesteps', or 'BD'
timeaxis = 'BD'

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = 
df = pd.read_csv(data_directory+'/networkx-allframes.csv')

ncolloids = pd.unique(df['ncolloids'])
frames = df['frame'].to_numpy()
lcc = df['lcc_size'].to_numpy() / ncolloids * 100

# set labels:
if timeaxis=='frame':
  times = frames
  max_time = max(times)
  axis_text = "Number of frames"
if timeaxis=='timesteps':
  times = frames*period
  max_time = max(times)
  axis_text = 'Number of timesteps'
if timeaxis=='BD':
  t1 = period * dt_Integration
  times = frames*t1
  max_time = max(times)
  axis_text = 'Brownian Dynamics (BD) times'

# plot data
plt.plot(times, lcc, 'grey', label="colloids") 

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.ylim([0, 105]) # max 100%
plt.xlim([0, max_time]) 
#plt.xticks(np.arange(0, max_time+1, 2.0))

plt.title('Largest Connected Component', fontsize=20)

plt.xlabel(axis_text, fontsize=16)          
plt.ylabel('Percent of colloids in LCC', fontsize=16)

plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'% at $D_0$='+str(D0)+'kT', title_fontsize=12) #loc=7, bbox_to_anchor=(1,0.85))
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

plt.savefig('lcc_phi'+str(phi)+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('LCC plot created')
