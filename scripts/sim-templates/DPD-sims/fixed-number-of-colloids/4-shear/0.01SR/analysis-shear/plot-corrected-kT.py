# plot the (corrected) kinetic temperature (kT) of a simulation
# (from GSD file) versus DPD time
# NOTE: requires "T_corrected.txt" created with the
#	sim-analysis-shear.py script

import numpy as np
import matplotlib.pyplot as plt

kT = 0.1

f = np.genfromtxt('T_corrected.txt')

# number of timesteps: frame*period
period = 10000
timesteps = f[:,0]*period

# DPD time: frame*t1
dt_Integration = 0.001
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1

plt.plot(f[:,0],f[0:,4])
plt.axhline(y=kT, color="r")

plt.xlabel('DPD times')
plt.ylabel('Corrected kT')

plt.yscale('log')
plt.xscale('log')

#plt.savefig('corrected-kT.png',dpi=600, transparent=False)
plt.show()
