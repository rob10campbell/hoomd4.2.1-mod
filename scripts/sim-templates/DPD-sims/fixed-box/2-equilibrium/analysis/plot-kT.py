# plot the kinetic temperature (kT) of a simulation
# (from GSD file) versus DPD time
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis.py script

import numpy as np
import matplotlib.pyplot as plt

kT = 0.1

f = np.genfromtxt('gsd-properties.txt')

# number of timesteps: frame*period
period = 10000
timesteps = f[:,0]*period

# DPD time: frame*t1
dt_Integration = 0.001
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1

plt.plot(DPDtimes,f[0:,8])
plt.axhline(y=kT, color="r")

plt.xlabel('DPD times')
plt.ylabel('kT')

plt.yscale('log')
plt.xscale('log')

#plt.savefig('kT.png',dpi=600, transparent=False)
plt.show()
