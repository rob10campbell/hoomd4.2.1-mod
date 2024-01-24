# plot the void/pore size distribution of a colloidal gel (from GSD file)
# options for both Torquato's Pore Size Distribution and Gubbin's 
# Pore Size Distribution versus void radius (r)
# NOTE: requires "poresize.txt" created with the sim-analysis-DPD.py, 
#       module_analysis_DPD.f90, and solvopt.f90 scripts

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

# choose timeaxis for plot:  'frame', 'timesteps', or 'DPD'
timeaxis = 'DPD'

# select type of pore size to plot (default: Gubbin's PSD)
#sizetype = 'Torquato' # volume where you can place an entire probe particle
sizetype = 'Gubbin' # volume occupied by any part of a probe particle

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = frame, probe_posx, probe_posy, probe_posz, 
#     voidcenter_x, voidcenter_y, voidcenter_z, porediameter_T, porediameter_G
pore_df_allframes = pd.read_csv(data_directory+'/poresize.csv')

# get the number of frames in the data (starts at 0)
nframes = len(pd.unique(pore_df_allframes['frame']))

# convert negative values into specific frames 
for i in range(len(framechoice)):
  if framechoice[i] < 0:
    framechoice[i] = nframes - abs(framechoice[i])

# select poresize only for the desired frames 
pore_df_framechoice = pore_df_allframes.loc[pore_df_allframes['frame'].isin(framechoice)][['frame','porediameter_T','porediameter_G']]

# find the maximum pore-size bin in the data
max_T = max(pore_df_framechoice['porediameter_T'])
max_G = max(pore_df_framechoice['porediameter_G'])
largest_pore = int(np.ceil(max(max_T, max_G)))

# get the number of pores in the last frame
lastframe_npores = len(pore_df_framechoice.loc[pore_df_framechoice['frame'] == framechoice[-1]]['porediameter_T'])

# get colors
n_curves = len(framechoice)
values = range(n_curves)
jet = cm = plt.get_cmap('plasma')
cNorm  = colors.Normalize(vmin=1, vmax=values[-1]+1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# add phi and D0 label to legend
plt.plot([], [], ' ', label='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT')


# bin and plot data for chosen frames
for i in range(len(framechoice)):

  # select data for the current frame
  frame = framechoice[i]
  pore_df_frame = pore_df_framechoice.loc[pore_df_framechoice['frame'] == frame]

  pore_diameters_T = pore_df_frame['porediameter_T'].to_numpy()
  pore_diameters_G = pore_df_frame['porediameter_G'].to_numpy()

  npores = len(pore_diameters_T)
  if npores != lastframe_npores:
    print('Warning: frame '+str(framechoice[i])+' and frame '+str(framechoice[-1])+' do not have the same number of probe points')

  # create poresize bins
  diameter_bins = np.linspace(0,largest_pore,int(largest_pore+1)).astype(int)

  # make array to count voids per bin
  histbins_T = np.zeros(int(largest_pore+1))
  histbins_G = np.zeros(int(largest_pore+1))

  # make array for data normalized by total number of voids
  normcounts_T = np.zeros(int(largest_pore+1))
  normcounts_G = np.zeros(int(largest_pore+1))

  # bin the voids
  for curr_bin in range(0,int(largest_pore)):
    for j in range(npores):
      if curr_bin <= pore_diameters_T[j] < curr_bin+1:
         histbins_T[curr_bin] += 1
    for j in range(npores):
      if curr_bin <= pore_diameters_G[j] < curr_bin+1:
         histbins_G[curr_bin] += 1

  normcounts_T[:] = histbins_T[:]/npores
  normcounts_G[:] = histbins_G[:]/npores

  # convert frames to times
  if timeaxis=='timesteps':
    # to plot vs. number of timesteps: frame*period
    time = frame*period
  if timeaxis=='DPD':
  # to plot vs. DPD time: frame*t1
    t1 = period * dt_Integration
    time = frame*t1

  # set colors:
  symcolor = scalarMap.to_rgba(values[i])
  # or match colors from ref paper (https://doi.org/10.1021/acs.macromol.9b02166)
  # Torquato's in pale blue: #a7b1ff
  # Gubbin's in pale red: #f6aea7

  # set labels:
  if timeaxis=='frame':
    labeltext = "frame "+str(frame)
  if timeaxis=='timesteps':
    labeltext = str(int(time))+' timesteps'
  if timeaxis=='DPD':
    labeltext = str(int(time))+' DPD times'

  if sizetype=='Torquato':
    plt.plot(diameter_bins/2, normcounts_T, color=symcolor, linestyle='--', label=labeltext)
    #plt.plot(diameter_bins/2, histbins_T, color='#a7b1ff', linestyle='-', label=labeltext)
    #plt.bar(diameter_bins/2, histbins_T, color='#a7b1ff', label=labeltext)
    legend_title="Torquato's PSD"

  if sizetype=='Gubbin':
    plt.plot(diameter_bins/2, normcounts_G, color=symcolor, linestyle='-', label=labeltext)
    #plt.plot(diameter_bins/2, histbins_G, color='#f6aea7', linestyle='-', label=labeltext)
    #plt.bar(diameter_bins/2, histbins_G, color='#f6aea7', label=labeltext)
    legend_title="Gubbin's PSD"

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

#plt.yscale('log')
#plt.xscale('log')

plt.title('Void Size Distribution ('+str(lastframe_npores)+' probe points)', fontsize=16)
plt.xlabel('Void Radius $r$ [DPD units]', fontsize=16) 
plt.ylabel('$P(d)$', fontsize=16)
#plt.ylabel('Number of Voids', fontsize=16)

plt.legend(prop={"size":12}, title=legend_title, title_fontsize=12) #,loc=7, bbox_to_anchor=(1,0.4))
plt.rcParams['figure.figsize'] = [6, 6]
plt.tight_layout()

plt.savefig('void-dist_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('pore size (void size) distribution plot created')
