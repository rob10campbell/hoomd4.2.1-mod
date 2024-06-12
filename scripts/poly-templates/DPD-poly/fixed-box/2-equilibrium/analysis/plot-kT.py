# plot the kinetic temperature (kT) of a simulation
# (from GSD file) versus DPD time
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis.py script

import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt('gsd-properties.txt')
plt.plot(f[:,0],f[0:,8])
plt.xlabel('DPD times')
plt.ylabel('kT')

kT = 0.1
plt.axhline(y=kT, color="r")

plt.yscale('log')
plt.xscale('log')

#plt.savefig('kT.png',dpi=600, transparent=False)
plt.show()
