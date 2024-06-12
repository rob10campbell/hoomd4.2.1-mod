# plot the contact number distribution for a single frame 
# of a BD colloidal sim
#    probability vs. coordination number
# NOTE: requires "Z-counts.csv" created with the
#	sim-analysis-BD.py and module_analysis_BD.f90 scripts

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

# choose which frames to plot (default last frame)
framechoice = [-1]
#framechoice = [0,10,20,50,100,150]

# choose timeaxis for plot:  'frame', 'timesteps', or 'BD'
timeaxis = 'BD'

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = frame, colloidID, typeid, Z 
Zdf = pd.read_csv(data_directory+'/Z-counts.csv')

# FOR ALL FRAMES find the unique frames so we can loop through them
#frames = pd.unique(Zdf['frame'])
#nframes = len(frames)

# FOR FRAMECHOICE
# convert negative values into specific frames 
nframes = len(pd.unique(Zdf['frame']))
for i in range(len(framechoice)):
  if framechoice[i] < 0:
    framechoice[i] = nframes - abs(framechoice[i])
frames = framechoice
nframes = len(framechoice)

# get the total number of colloids from the last frame
last_frame = Zdf.loc[Zdf['frame'] == frames[-1]]
ncolloids = max(last_frame['colloidID'])

# make a list of bins for the x-axis
max_Z = max(Zdf['Z'].to_numpy())
bins = np.linspace(0,max_Z,int(max_Z)).astype(int)

# create empty bins for the Z values
histbins = np.zeros((nframes,max_Z))
probs = np.zeros((nframes,max_Z))

# get colors
n_curves = len(framechoice)
values = range(n_curves)
jet = cm = plt.get_cmap('plasma')
cNorm  = colors.Normalize(vmin=1, vmax=values[-1]+1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


for f in range(nframes):
  frame = frames[f]
  # collect the contact numbers for this frame in each group
  Zs_frame = Zdf.loc[Zdf['frame'] == frame]['Z'].to_numpy() 

  # bin the Zs
  for i in range(1,max_Z+1):
    for z in Zs_frame:
       if z == i:
         histbins[f][i-1] += 1
 
  # convert to probability
  probs[f] = histbins[f]/ncolloids 

  # set labels:
  if timeaxis=='frame': 
    labeltext = "frame "+str(frame)
  if timeaxis=='timesteps':
    time = frame*period
    labeltext = str(int(time))+' timesteps'
  if timeaxis=='BD':
    t1 = period * dt_Integration
    time = frame*t1
    labeltext = str(int(time))+' BD times'

  symcolor = scalarMap.to_rgba(values[f])    
  plt.plot(bins, probs[f], color=symcolor, label=labeltext)
  #plt.bar(bins, probs[f], color=symcolor, alpha=1, label=labeltext)

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.title('Contact Number Distribution', fontsize=16)
plt.xlabel('$Z$', fontsize=16)
plt.ylabel('$P(Z)$', fontsize=16)

#ax.set_xticks(range(0,max_Z,2))

# figsize includes title, axes, and plot; move lgd before this to include lgd in figsize
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

lgd = plt.legend(prop={"size":12}, title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12, loc='center left', bbox_to_anchor=(1,0.5))

plt.savefig('Zdist_phi'+str(int(phi))+'_'+str(D0)+'kT.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600, transparent=False)
#plt.show() # NOTE: plot.show() will cut off the legend, but the figure saves correctly

print('Z distribution plot created')
 
