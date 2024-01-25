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
fig.suptitle('Radial Distribution Function g(r) ($\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT)', fontsize=16)

axs[0].axhline(y=1, linewidth=1, linestyle='--', color='black') 
axs[0].plot(f[:,0], f[:,1], 'royalblue', label='colloid-solvent')
#axs[0].set_title('colloid-solvent', fontsize=16)
# add particle size to the legend
axs[0].plot([], [], ' ', label='$r_{C} = $'+str(R_C))


axs[1].axhline(y=1, linewidth=1, linestyle='--', color='black') 
axs[1].plot(f[:,0], f[:,2], 'dimgrey', label='colloid-colloid')
#axs[1].set_title('colloid-colloid', fontsize=16)
# add particle size to the legend
axs[1].plot([], [], ' ', label='$r_{C} = $'+str(R_C))

for ax in axs.flat:
    ax.set(xlabel='$r$ [DPD units]', ylabel='$g(r)$')
    ax.legend(prop={"size":12})
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)

plt.savefig('gofr_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('g(r) plot created')
