# plot the radial distribution function g(r) for colloids in
# the final frame of a DPD sim (from GSD file)
#    this the number of particles per unit volume versus 
#    DPD distance (distance normalized by R_C)
# NOTE: requires "gofr.txt" created with the
#	sim-analysis-DPD.py and module_analysis_DPD.f90 scripts

import numpy as np
import matplotlib.pyplot as plt

##############
""" INPUTS """
##############

# sim parameters
phi = 20
D0 = 12
R_C = 1 

# data source
data_directory = '../data'

##########################
""" ANALYSIS AND PLOT  """
##########################

# load data: columns = position gofr_cs gofr_cc
f = np.genfromtxt(data_directory+'/gofr.txt', skip_header=1)

fig, axs = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

axs[0].axhline(y=1, linewidth=1, linestyle='--', color='black') 
axs[0].plot(f[:,0], f[:,1], 'royalblue', label='colloid-solvent')
axs[0].plot([], [], ' ', label='colloid-solvent')
#axs[0].set_title('colloid-solvent', fontsize=16)

# label subplot with custom sub-legend 
axs[0].plot([], [], ' ', label='colloid-solvent')
axs[0].plot([], [], ' ', label='$r_{C} = $'+str(R_C))
lines_cs = axs[0].get_lines()
lgd_labels_cs = axs[0].legend(prop={"size":12}, handles=[lines_cs[i] for i in [-2,-1]], labels=[lines_cs[i].get_label() for i in [-2,-1]], loc='upper right')
axs[0].add_artist(lgd_labels_cs)


axs[1].axhline(y=1, linewidth=1, linestyle='--', color='black') 
axs[1].plot(f[:,0], f[:,2], 'dimgrey', label='colloid-colloid')
#axs[1].set_title('colloid-colloid', fontsize=16)

# label subplot with custom sub-legend 
axs[1].plot([], [], ' ', label='colloid-colloid')
axs[1].plot([], [], ' ', label='$r_{C} = $'+str(R_C))
lines_cc = axs[1].get_lines()
lgd_labels_cc = axs[1].legend(prop={"size":12}, handles=[lines_cc[i] for i in [-2,-1]], labels=[lines_cc[i].get_label() for i in [-2,-1]], loc='upper right')
axs[1].add_artist(lgd_labels_cc)

for ax in axs.flat:
    ax.set(xlabel='$r$ [DPD units]', ylabel='$g(r)$')
    #ax.legend(prop={"size":12})
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)

# add one legend for the lines in both subplots
lines_labels = ([lines_cs[i] for i in range(1,nseeds+2)],[lines_cs[i].get_label() for i in range(1,nseeds+2)]) # data only (same on both plots, so use axs[0])
lgd = fig.legend(prop={"size":12}, handles=lines_labels[0], labels=lines_labels[1], title='$\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT', title_fontsize=12, loc='center left', bbox_to_anchor=(1,0.5))

title = fig.text(0.5, 1.03, 'Radial Distribution Function g(r)', ha='center', va='center', fontsize=16)

plt.savefig('gofr_phi'+str(int(phi))+'_'+str(D0)+'kT.png',  bbox_extra_artists=(lgd,title), bbox_inches='tight', dpi=600, transparent=False)
#plt.show()

print('g(r) plot created')
