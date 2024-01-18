# plot the distribution of voronoi volumes for a single frame
# in a biomdal sim (from GSD file)
#    probability vs. voronoi volume
# NOTE: requires "voronoi_binned.csv" created with the
#	sim-voronoi-DPD.py script

import numpy as np
import pandas as pd
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

##############
""" INPUTS """
##############

# sim parameters 
phi = 20
D0 = 12
R_C = 1
period = 10000
dt_Integration = 0.001
sim_vol = 70*70*70 

# choose timeaxis for plot:  'frame', 'timesteps', or 'DPD'
timeaxis = 'DPD'

# choose which frames to plot (default last frame)
#framechoice = [-1]
framechoice = [0,10,20,50,100,150]

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = vol-bins, mono-counts, mono-norm
filepath = data_directory+'/voronoi_binned.csv'
vorodist_df = pd.read_csv(filepath)

# get the number of frames in the data (starts at 0)
frames = pd.unique(vorodist_df['frame'])
nframes = len(pd.unique(vorodist_df['frame']))

# convert negative values into specific frames 
for i in range(len(framechoice)):
  if framechoice[i] < 0:
    framechoice[i] = nframes - abs(framechoice[i])

# create bins to store counts at each volume
bins = np.unique(vorodist_df['vol-bins'].to_numpy())

# create a color gradient to represent time
def create_color_gradient(start_color, end_color, num_colors):
  cmap = mcolors.LinearSegmentedColormap.from_list(
      name='custom_cmap',
      colors=[start_color, end_color],
      N=num_colors
  )
  return [cmap(i) for i in np.linspace(0, 1, num_colors)]
start_color = 'lightgrey'
end_color = 'black'
num_colors = len(framechoice)+2
# create gradient list
color_gradient = create_color_gradient(start_color, end_color, num_colors)


plt.title('Voronoi Volume Distribution', fontsize=20) #, pad=20)

for i in range(len(framechoice)):
  frame = framechoice[i]
  if frame > nframes:
    break
  frame_df = vorodist_df.loc[vorodist_df['frame'] == frame]

  # set time label:
  if timeaxis=='frame':
    time = frame
    label_text = "frame "+str(time)
  if timeaxis=='timesteps':
    time = frame
    label_text = str(time)+' timesteps'
  if timeaxis=='DPD':
    t1 = period * dt_Integration
    time = frame
    label_text = str(time)+' DPD times'

  plt.plot(bins, frame_df['mono-norm'], linestyle='-', color=color_gradient[i+1], label=label_text)

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.ylabel('$P(V_{voro})$', fontsize=16)
plt.xlabel('$V_{voro}$ [cubic DPD units]', fontsize=16)
#plt.xticks(range(0,true_max_Z+1,1))
plt.xlim(0.1,sim_vol)
plt.ylim(0.0001,1)

lgd = plt.legend(prop={"size":12}, loc="upper right", title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)
plt.tight_layout() 
plt.rcParams['figure.figsize'] = [6, 6]

plt.loglog()
plt.savefig('voronoi_loglog_phi'+str(int(phi))+'_'+str(D0)+'kT.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=False)

#plt.show()
plt.close()

print('Voronoi distribution plot created')
