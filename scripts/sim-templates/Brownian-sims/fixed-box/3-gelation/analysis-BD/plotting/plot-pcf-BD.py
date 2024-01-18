# plot the qualitative approximation of the pair correlation function 
# h(r) for colloids in a BD sim (from GSD file)
#    this is the probability density per unit volume of a 
#    particle at r1 finding another particle at r2 
#    (plotted as 3 2D planes)
# NOTE: requires "pcf.txt" created with the
#	sim-analysis-BD.py and module_analysis_BD.f90 scripts

import numpy as np
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, cm
from matplotlib import colors as CM
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

##############
""" INPUTS """
##############

# set data-specific parameters
phi = 20
D0 = 12
dims_max = 10
R_C = 1

# data source
data_directory = '../data'


##########################
""" ANALYSIS AND PLOT  """
##########################

# choose color scheme
colors = ['bwr', 'plasma', 'viridis', 'inferno', 'magma', 'cividis', 'Spectral', 'twilight']
color = colors[1]


# read in particle count data
df_counts_xy_cc = pd.read_csv(data_directory+'/PCF_counts_XY.csv')
counts_XY_CC = df_counts_xy_cc.to_numpy()
df_counts_xz_cc = pd.read_csv(data_directory+'/PCF_counts_XZ.csv')
counts_XZ_CC = df_counts_xz_cc.to_numpy()
df_counts_yz_cc = pd.read_csv(data_directory+'/PCF_counts_YZ.csv')
counts_YZ_CC = df_counts_yz_cc.to_numpy()

# set 2D plot dimensions
extent_set = [-dims_max, dims_max, -dims_max, dims_max]


### PLOT COUNTS (1 FIGURE WITH 3 PLOTS)
fig, (ax_xy, ax_xz, ax_yz) = plt.subplots(1,3, figsize=(9,3))
fig.suptitle('Pair Correlation Function ($\phi$='+str(phi)+'% at $D_0$='+str(D0)+'kT)', fontsize=16)
plt.setp( (ax_xy, ax_xz, ax_yz), 
  xticks=[-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], 
  yticks=[-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)] )

## XY Subplot
pcf_xy = ax_xy.imshow(counts_XY_CC,extent = extent_set, cmap=color)
ax_xy.set_title("XY Plane", fontsize="14")
ax_xy.set_xlabel('Separation Distance [$R_C$]', fontsize = 10)
ax_xy.set_ylabel('Separation Distance [$R_C$]', fontsize = 10)

ax_xy.set_xticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=8)
ax_xy.set_yticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=8)
#ax_xy.grid(None)

# show a colloid particle centered in the plot
#circle_xy = plt.Circle((0, 0), R_C, color='white', alpha=0.5)
#ax_xy.add_patch(circle_xy)

# add colorbar
divider_xy = make_axes_locatable(ax_xy)
cax_xy = divider_xy.append_axes("right", size="5%", pad=0.1) 
cbar_xy = fig.colorbar(pcf_xy, cax=cax_xy)
cbar_xy.ax.text(3.5, 1.05, '# colloids', ha='center', va='center', fontsize=10, transform=cbar_xy.ax.transAxes)
cbar_xy.locator = MaxNLocator(nbins=5)
cbar_xy.update_ticks()
cbar_xy.ax.tick_params(axis='y', labelsize=10)


## XZ Subplot
pcf_xz = ax_xz.imshow(counts_XZ_CC,extent = extent_set, cmap=color)
ax_xz.set_title("XZ Plane", fontsize=14)
ax_xz.set_xlabel('Separation Distance [$R_C$]', fontsize = 10)
#ax_xz.set_ylabel('Separation Distance [$R_C$]', fontsize = 10)
ax_xz.set_xticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=8)
ax_xz.set_yticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=8)
#ax2_xz.grid(None)

# show a colloid particle centered in the plot
#circle_xz = plt.Circle((0, 0), 1, color='white', alpha=0.5)
#ax_xz.add_patch(circle_xz)

# add colorbar
divider_xz = make_axes_locatable(ax_xz)
cax_xz = divider_xz.append_axes("right", size="5%", pad=0.1) 
cbar_xz = fig.colorbar(pcf_xz, cax=cax_xz)
cbar_xz.ax.text(3.5, 1.05, '# colloids', ha='center', va='center', fontsize=10, transform=cbar_xz.ax.transAxes)
cbar_xz.locator = MaxNLocator(nbins=5) 
cbar_xz.update_ticks()
cbar_xz.ax.tick_params(axis='y', labelsize=10)


## YZ Subplot
pcf_yz = ax_yz.imshow(counts_YZ_CC,extent = extent_set, cmap=color)
ax_yz.set_title("YZ Plane", fontsize=14)
ax_yz.set_xlabel('Separation Distance [$R_C$]', fontsize = 10)
#ax_yz.set_ylabel('Separation Distance [$R_C$]', fontsize = 10)
ax_yz.set_xticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=10)
ax_yz.set_yticklabels([-int(dims_max), -int(dims_max/2), 0, int(dims_max/2), int(dims_max)], fontsize=10)
#ax_yz.grid(None)

# show a colloid particle centered in the plot
#circle_yz = plt.Circle((0, 0), 1, color='white', alpha=0.5)
#ax3.add_patch(circle_yz)

# add colorbar
divider_yz = make_axes_locatable(ax_yz)
cax_yz = divider_yz.append_axes("right", size="5%", pad=0.1) 
cbar_yz = fig.colorbar(pcf_yz, cax=cax_yz)
cbar_yz.ax.text(3.5, 1.05, '# colloids', ha='center', va='center', fontsize=10, transform=cbar_yz.ax.transAxes)
cbar_yz.locator = MaxNLocator(nbins=5)
cbar_yz.update_ticks()
cbar_yz.ax.tick_params(axis='y', labelsize=10)

# save/show figure
fig.tight_layout()
plt.savefig('PCF_counts_phi'+str(int(phi))+'_'+str(D0)+'kT.png', dpi=1000, transparent=False)
#plt.show()
plt.close()

print('PCF plot created')
