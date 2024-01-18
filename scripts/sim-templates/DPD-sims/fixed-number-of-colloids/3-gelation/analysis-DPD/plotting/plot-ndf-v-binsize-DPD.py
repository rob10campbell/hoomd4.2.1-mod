# plot the binned number density fluctuation of all colloids 
# in a DPD sim (from GSD file) for different bin sizes; there
# should be a peak at the average cluster size
#     Index of dispersion (D or VMR) vs. time
#     (time can be frames, timesteps, or DPD time)
# NOTE: requires "ndfluc.txt" created with the
#	sim-analysis-DPD.py script

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
period = 10000
dt_Integration = 0.001

# optional: plot from a limited size-range (see plotting section below)
#plot_sizerange_start = 7

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: colums = "frame" "bin-size" "ndfluc_C"
df = pd.read_table(data_directory+'/ndfluc.txt', delimiter=" ")

def plot_ndfluc_py(flag):

  # pick all frames OR one frame
  #for frame in pd.unique(df['frame']): # for all frames
    frame = max(pd.unique(df['frame'])) # for one frame

    frame_df = df.loc[df['frame']==frame]

    # to plot vs. simulation frame number
    if flag=='frame':
      time = frame 
   
    # to plot vs. number of timesteps: frame*period
    if flag=='step':
      time = frame*period

    # to plot vs. DPD time: frame*t1
    if flag=='DPD':
      t1 = period * dt_Integration
      time = frame*t1


    # create arrays to store data for this frame
    bins = pd.unique(df['bin-size'])
    nbins = len(bins)
    ndf_all_bybin = []

    # get the data for this frame    
    for bin_size in pd.unique(df['bin-size']):
      ndf_all = frame_df.loc[frame_df['bin-size']==bin_size]['ndfluc_C'].to_numpy()[0]
      ndf_all_bybin.append(ndf_all)
 
    # create the plot for this frame
    plt.rcParams['figure.figsize'] = [6, 6]

    # select lable to match the units of time 
    if flag=='frame':
      label_text = 'frame '+str(round(time))
    if flag=='step':
      label_text = str(round(time))+' timesteps'
    if flag=='DPD':
      label_text = str(round(time))+' DPD times'

    plt.plot(bins,ndf_all_bybin, color='grey', marker='o', linestyle='-', label=label_text)
    # OR, optional: plot from a limited size-range with no markers
    #plt.plot(bins[plot_sizerange_start-1:],ndf_all_bybin[plot_sizerange_start-1:], 'grey', linestyle='-', label=label_text)

    # add number of bins to the legend
    plt.plot([], [], ' ', label=str(nbins)+' bin sizes')

    #plt.yscale('log')
    #plt.xscale('log')
    #plt.ylim(0,1)
    plt.xticks(np.arange(0, max(bins)+1, 5))

    plt.title('Number Density Fluctuation by Bin Size', fontsize=16)
    plt.xlabel('bin size [DPD units]')
    plt.ylabel('$(\langle N^2 \\rangle - \langle N \\rangle^2) / \langle N \\rangle$', fontsize=16)

    
    plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)# loc='upper right')

    plt.tight_layout()

    plt.savefig('ndfluc-vs-binsize-frame'+str(round(frame))+'_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False) 
    #plt.show()
    plt.close()

    print('NDF vs. bin size plot created')

if __name__ == '__main__':
  #plot_ndfluc_py('frame')
  #plot_ndfluc_py('step')
  plot_ndfluc_py('DPD')
