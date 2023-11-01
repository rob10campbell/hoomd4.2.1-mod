# semi-log plot of (corrected) stress versus strain
# NOTE: requires "T_corrected.txt" created with the
#	sim-analysis-shear.py script

import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt('T_corrected.txt')

# number of timesteps: frame*period
period = 10000
timesteps = f[:,0]*period

# DPD time: frame*t1
dt_Integration = 0.001
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1

plt.plot(DPDtimes,-f[0:,1]) #xx
#plt.plot(DPDtimes,-f[0:,2]) #xy
#plt.plot(DPDtimes,-f[0:,3]) #xz

plt.xscale('log')

plt.xlabel('DPD times')
plt.ylabel('Corrected Shear Stress')

#plt.savefig('corrected-stress.png',dpi=600, transparent=False)
plt.show()
