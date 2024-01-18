# plot the kurtosis and/or excess kurtosis (4th moment) of 
# voronoi volumes for all frames of a BD colloid sim (from GSD)
# NOTE: requires "voronoi-stats.csv" created with the
#       sim-voronoi-BD.py script

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

# choose timeaxis for plot:  'frame', 'timesteps', or 'BD'
timeaxis = 'BD'

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = frame, mean, divergence, skew, kurtosis, excess-kurtosis
stats_df = pd.read_csv(data_directory+'/voronoi-stats.csv')

def plot_voronoi_kurtosis(stats_df, bare=False, excess=True, show=False):
  # default behavior:
  #  - plots Fisher excess kurtosis and not the Pearson bare kurtosis
  #  - saves plots (not show)

  # extract data for plotting
  frame_numbers = stats_df['frame'].to_numpy()
  kurtosis_volume = stats_df['kurtosis'].to_numpy()
  excesskurtosis_volume = stats_df['excess-kurtosis'].to_numpy()

  # set time labels:
  if timeaxis=='frame':
    time = frame_numbers
    axis_text = "Frame Number"
  if timeaxis=='timesteps':
    time = frame_numbers*period
    axis_text = 'Timesteps'
  if timeaxis=='BD':
    t1 = period * dt_Integration
    time = frame_numbers*t1
    axis_text = 'Brownian Dynamics (BD) times'

  # 4. Plot the outlier-ness of voronoi volumes
  # Kurtosis is a measure of the tails -- are there more outliers than in a normal distribution?
  # < 3 == thin tails; fewer and/or less extreme outliers than the normal distribution
  # > 3 == fat tails; more outliers than the normal distribution

  # a. Kurtosis (Pearson)
  if bare==True: 
    plt.plot(time, kurtosis_volume, linestyle='-', color='grey')

    # add particle size to the legend
    plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

    plt.title('Are There More Outliers Than "Normal"?', fontsize=16)
   
    plt.ylabel('Bare Kurtosis', fontsize=16)
    plt.xlabel(axis_text, fontsize=16)

    plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
    plt.rcParams['figure.figsize'] = [6, 6]
    plt.tight_layout()
    #plt.grid(True)

    if show==True:
      plt.show()
    else:
      plt.savefig('voronoi-kurtosis_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)

    print('Voronoi dist kurtosis plot created')

  # b. Excess Kurtosis (Fisher)
  # NEGATIVE == thin tails; fewer and/or less extreme outliers than the normal distribution
  # POSITIVE == fat tails; more outliers than the normal distribution
  if excess==True: 
    plt.plot(time, excesskurtosis_volume, linestyle='-', color='grey')

    # add particle size to the legend
    plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

    plt.title('Are There More Outliers Than "Normal"?', fontsize=16)
    
    plt.ylabel('Excess Kurtosis', fontsize=16)
    plt.xlabel(axis_text, fontsize=16)

    plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
    plt.rcParams['figure.figsize'] = [6, 6]
    plt.tight_layout()
    #plt.grid(True)

    if show==True:
      plt.show()
    else:
      plt.savefig('voronoi-excesskurtosis_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)

    print('Voronoi dist excess kurtosis plot created')

if __name__ == '__main__':
   plot_voronoi_kurtosis(stats_df)
   #plot_voronoi_kurtosis(stats_df, show=True)
   #plot_voronoi_kurtosis(stats_df, bare=True, excess=False, show=False)
   #plot_voronoi_kurtosis(stats_df, bare=True, excess=False, show=True) 
