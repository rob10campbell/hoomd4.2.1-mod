# plot the number density fluctuation of colloids and solvent
# in a BD sim (from GSD file)
#     Index of dispersion (D or VMR) vs. time
#     (time can be frames, timesteps, or BD time)
# NOTE: requires "ndfluc.txt" created with the
#	sim-analysis-BD.py script

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

##############
""" INPUTS """
##############

# sim parameters
phi = 20
D0 = 12
R_C = 1
period = 10000
dt_Integration = 0.001

# choose which bin sizes to plot (default size 10)
binchoice = [10]
#binchoice = [2,5,10,15,35,70]

# choose timeaxis for plot:  'frame', 'timesteps', or 'BD'
timeaxis = 'BD'

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: column names = frame bin-size ndfluc_C
df = pd.read_table(data_directory+'/ndfluc.txt', delimiter=" ")
frames = pd.unique(df['frame'])
bin_sizes = pd.unique(df['bin-size'])

# limit data to the chosen bins
df_binchoice = df.loc[df['bin-size'].isin(binchoice)]
nbinsizes = len(binchoice)

# convert frames to times
if timeaxis=='frame':
  time = frames
  axis_text = "Number of frames"
if timeaxis=='timesteps':
  time = frames*period
  axis_text = 'Number of timesteps'
if timeaxis=='BD':
  t1 = period * dt_Integration
  time = frames*t1
  axis_text = 'Brownian Dynamics (BD) times'

# get colors
n_curves = nbinsizes
values = range(n_curves)
jet = cm = plt.get_cmap('plasma')
cNorm  = colors.Normalize(vmin=1, vmax=values[-1]+1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))
                                 
for i in range(nbinsizes):
  curr_bin = binchoice[i]

  bin_df = df.loc[df['bin-size'] == curr_bin]

  symcolor = scalarMap.to_rgba(values[i])
  plt.plot(time,bin_df['ndfluc_C'], color=symcolor, label="bin-size "+str(curr_bin)) 

plt.yscale('log')
plt.xscale('log')

plt.title('Number Density Fluctuation', fontsize=16)

plt.xlabel(axis_text, fontsize=16)
plt.ylabel('$\langle N^2 \\rangle - \langle N \\rangle^2 / \langle N \\rangle$', fontsize=16)

plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12) #, loc=7, bbox_to_anchor=(1,0.4))
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

plt.savefig('ndfluc_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('NDF plot created')
