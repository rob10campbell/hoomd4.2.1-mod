# plot the largest connected compoinent (LCC) of colloids 
# in DPD sim (from networkX analysis)
#    percent of colloids in the largest connected component versus time
#    (time can be frames, timesteps, or DPD time)
# NOTE: requires "networkx-allframes.csv" created with the
#	sim-networkCSV-DPD.py script

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

# choose timeaxis for plot:  'frame', 'timesteps', or 'DPD'
timeaxis = 'DPD'

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
if timeaxis=='DPD':
  t1 = period * dt_Integration
  times = frames*t1
  max_time = max(times)
  axis_text = 'DPD times'

# plot data
plt.plot(times, lcc, 'grey', label="colloids") 

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

# linear axes
#plt.ylim([0, 105]) # max 100%
#plt.xlim([0, max_time]) 
#plt.xticks(np.arange(0, max_time+1, 2.0))

# log axes
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-1, 125]) # max 100%
plt.xlim([1e-1, max_time])
#plt.xticks(np.arange(0, max_time+1, 2.0))

plt.title('Largest Connected Component', fontsize=20)

plt.xlabel(axis_text, fontsize=16)          
plt.ylabel('Percent of colloids in LCC', fontsize=16)

# figsize includes title, axes, and plot; move lgd before this to include lgd in figsize
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

lgd = plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12, loc='center left', bbox_to_anchor=(1,0.5)) #loc=7, bbox_to_anchor=(1,0.85))

plt.savefig('lcc_phi'+str(phi)+'_'+str(D0)+'kT.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600, transparent=False)
#plt.show() # NOTE: plt.show() will cut off the legend, but figure will save correctly

print('LCC plot created')
