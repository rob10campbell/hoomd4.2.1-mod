# plot the mean and divergence (standard deviation) of voronoi volumes 
# for all frames of a DPD colloid sim (from GSD file)
# NOTE: requires "voronoi-stats.csv" created with the
#       sim-voronoi-DPD.py script

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import skew
from scipy.stats import kurtosis
import sys

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

# load data: columns = frame, mean, divergence, skew, kurtosis, excess-kurtosis
stats_df = pd.read_csv(data_directory+'/voronoi-stats.csv')

# extract data for plotting
frame_numbers = stats_df['frame'].to_numpy()
mean_volume = stats_df['mean'].to_numpy()
sd_volume = stats_df['divergence'].to_numpy()

# set time labels:
if timeaxis=='frame':
  time = frame_numbers
  axis_text = "Frame Number"
if timeaxis=='timesteps':
  time = frame_numbers*period
  axis_text = 'Timesteps'
if timeaxis=='DPD':
  t1 = period * dt_Integration
  time = frame_numbers*t1
  axis_text = 'DPD times'


# 1. Plot the mean volume
# mean should be constant (constant total volume, constant total number of particles)

#'''
plt.plot(time, mean_volume, linestyle='-', color='grey')
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.title('Mean Voronoi Volume Over Time', fontsize=16)

plt.xlabel(axis_text, fontsize=16)
plt.ylabel('Mean Voronoi Volume', fontsize=16)

plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()
#plt.grid(True)

plt.savefig('voronoi-mean_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()
plt.close()

print('Voronoi dist mean plot created')
#'''

# 2. Plot the divergence of the volume distribution
plt.plot(time, sd_volume, linestyle='-', color='grey')

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.title('Divergence of Voronoi Volumes Over Time', fontsize=16)

plt.xlabel(axis_text, fontsize=16)
plt.ylabel('Divergence of Voronoi Volume', fontsize=16)

plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()
#plt.grid(True)

plt.savefig('voronoi-sd_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('Voronoi dist standard deviation plot created')

