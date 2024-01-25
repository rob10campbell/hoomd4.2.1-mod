# plot the radial distribution function g(r) for colloids in
# the final frame of a BD sim (from GSD file)
#    this the number of particles per unit volume versus 
#    BD distance (distance normalized by R_C)
# NOTE: requires "gofr.txt" created with the
#	sim-analysis-BD.py and module_analysis_BD.f90 scripts

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

# load data: columns = position gofr_cc
f = np.genfromtxt(data_directory+'/gofr.txt', skip_header=1)

plt.title('Radial Distribution Function g(r) ($\phi$='+str(phi)+'%, $D_0$='+str(D0)+'kT)', fontsize=16)

plt.axhline(y=1, linewidth=1, linestyle='--', color='black') 
plt.plot(f[:,0], f[:,1], 'dimgrey', label='colloid-colloid')

# add particle size to the legend
plt.plot([], [], ' ', label='$r_{C} = $'+str(R_C))

plt.xlabel('$r$ [BD units]', fontsize=16)
plt.ylabel('$g(r)$', fontsize=16)
plt.legend(prop={"size":12})

plt.savefig('gofr_phi'+str(int(phi))+'_'+str(D0)+'kT.png',dpi=600, transparent=False)
#plt.show()

print('g(r) plot created')
