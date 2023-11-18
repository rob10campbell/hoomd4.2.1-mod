# semi-log plot of (corrected) stress versus strain
# NOTE: requires "T_corrected.txt" created with the
#	sim-analysis-shear.py script

import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt('T_corrected.txt')


# NOTE: for shear analysis, we want time in terms of strain
# 	i.e. the trigger is "nframe_strain" not "period"
SR = 0.1 # from simulation
N_strains = 10 # from simulation
theta = 1.0 # from simulation
dt_Integration = 0.001 # from simulation
delta_T_shearing = round(theta/SR/dt_Integration)
nframe_strain = delta_T_shearing/10 # 20 frames per for each theta strain(s)
t1 = theta / nframe_strain * N_strains # timestep conversion to strains
strains = f[:,0]*t1


plt.plot(strains,-f[0:,1]) #xx
#plt.plot(strains,-f[0:,2]) #xy
#plt.plot(strains,-f[0:,3]) #xz

plt.xscale('log')

plt.xlabel('Strains')
plt.ylabel('Corrected Shear Stress')

#plt.savefig('corrected-stress.png',dpi=600, transparent=False)
plt.show()
