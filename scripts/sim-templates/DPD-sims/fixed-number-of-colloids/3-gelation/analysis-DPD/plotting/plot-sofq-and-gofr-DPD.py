# plot S(q) for colloids in a DPD sim (from GSD file)
# and optionally, plot the associated g(r)
#   S(q) is the structure factor calculated from the 
#     radial correlation function, h(r)=g(r)-1
#   g(r) is the radial distribution function, the number 
#     of particles per unit volume versus DPD distance 
#     (the distance normalized by R_C) 
# NOTE: requires "sofq.csv" and "gofr_sq.csv" created with the
#	sim-analysis-DPD.py and module_analysis_DPD.f90 scripts

import numpy as np
import pandas as pd
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

# choose timeaxis for plot:  'frame', 'timesteps', or 'DPD'
timeaxis = 'DPD'

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = 'frame' 'q' 'Sq' 
sq_df_allframes = pd.read_csv(data_directory+'/sofq.csv')

# get the number of frames in the data (starts at 0)
nframes = len(pd.unique(sq_df_allframes['frame']))

# get colors
n_curves = len(framechoice)
values = range(n_curves)
jet = cm = plt.get_cmap('plasma')
cNorm  = colors.Normalize(vmin=1, vmax=values[-1]+1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# optional: plot reference line at 1
#plt.axhline(y=1, linewidth=1, linestyle='--', color='black') 

# plot data for chosen frames
for i in range(len(framechoice)):

  if framechoice[i] < 0:
    frame = nframes - abs(framechoice[i])
  else:
    frame = framechoice[i]

  sq_df = sq_df_allframes[sq_df_allframes['frame'] == frame]

  # convert to desired time units
  if timeaxis=='frame':
    labeltext = "frame "+str(frame)
  if timeaxis=='timesteps': 
    time = frame*period
    labeltext = str(int(time))+' timesteps'
  if timeaxis=='DPD':
    t1 = period * dt_Integration
    time = frame*t1
    labeltext = str(int(time))+' DPD times'

  # plot the data
  symcolor = scalarMap.to_rgba(values[i])
  plt.plot(sq_df['q'], sq_df['Sq'], color=symcolor, label=labeltext)
  #plt.axhline(y=1, linewidth=1, linestyle='--', color='black') 

# add particle size to the legend
plt.plot([], [], ' ', label='$a = r_{C} = $'+str(R_C))

plt.title('Static Structure Factor $S(q)$', fontsize=16)
plt.xlabel('$qa$', fontsize=16)
plt.ylabel('$S(qa)$', fontsize=16)

plt.yscale('log')
plt.xscale('log')

#plt.ylim([0.4,15])
plt.xlim([0.08,2.2])

# figsize includes title, axes, and plot; move lgd before this to include lgd in figsize
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

lgd = plt.legend(prop={"size":10}, loc='center left', bbox_to_anchor=(1,0.5), title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)

plt.savefig('sofq_phi'+str(int(phi))+'_'+str(D0)+'kT.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600, transparent=False)
#plt.show() #NOTE: plot.show() will cut off the legend, but the figure will save correctly

plt.close()

print('S(q) plot created')

#'''
# plot matching g(r)

# load data: columns = 'frame', 'position', 'gofr' 
gofr_df_allframes = pd.read_csv(data_directory+'/gofr_sq.csv')

plt.axhline(y=1, linewidth=1, linestyle='--', color='black') 

for i in range(len(framechoice)):

  if framechoice[i] < 0:
    frame = nframes - abs(framechoice[i])
  else:
    frame = framechoice[i]

  gofr_df = gofr_df_allframes[gofr_df_allframes['frame'] == frame]

  # convert frames to times
  if timeaxis=='timesteps':
    # to plot vs. number of timesteps: frame*period
    time = frame*period
  if timeaxis=='DPD':
  # to plot vs. DPD time: frame*t1
    t1 = period * dt_Integration
    time = frame*t1

  # set labels:
  if timeaxis=='frame':
    labeltext = "frame "+str(frame)
  if timeaxis=='timesteps':
    labeltext = str(int(time))+' timesteps'
  if timeaxis=='DPD':
    labeltext = str(int(time))+' DPD times'

  # plot the data
  symcolor = scalarMap.to_rgba(values[i])
  plt.plot(gofr_df['position'], gofr_df['gofr'], color=symcolor, label=labeltext)

# add particle size to the legend
plt.plot([], [], ' ', label='$d_{C} = 2r_{C} = $'+str(2*R_C))

plt.title('Radial Distribution Function $g(r)$', fontsize=16)
plt.xlabel('$r$ [DPD units]', fontsize=16)
plt.ylabel('$g(r)$', fontsize=16)

# figsize includes title, axes, and plot; move lgd before this to include lgd in figsize
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

lgd = plt.legend(prop={"size":10}, loc='center left', bbox_to_anchor=(1,0.5), title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12)

plt.savefig('sq-gofr_phi'+str(int(phi))+'_'+str(D0)+'kT.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=600, transparent=False)
#plt.show() #NOTE: plot.show() will cut off the legend, but the figure will still save correctly

plt.close()

print('Plot of g(r) data used to calculate S(q) created')

#'''
