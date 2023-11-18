# plot the (corrected) kinetic temperature (kT) of a simulation
# (from GSD file) versus DPD time
# NOTE: requires "T_corrected.txt" created with the
#	sim-analysis-shear.py script

import numpy as np
import matplotlib.pyplot as plt

kT = 0.1

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


plt.plot(strains,f[0:,4])
plt.axhline(y=kT, color="r")

plt.xlabel('Strains')
plt.ylabel('Corrected kT')

plt.yscale('log')
plt.xscale('log')

#plt.savefig('corrected-kT.png',dpi=600, transparent=False)
plt.show()
