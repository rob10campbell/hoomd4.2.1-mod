# From a GSD file, plot the shear stress (-pressure_xy) 
# of a simulation versus DPD time
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis.py script

import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt('gsd-properties.txt', skip_header=1)

# number of timesteps: frame*period
period = 10000
timesteps = f[:,0]*period

# DPD time: frame*t1
dt_Integration = 0.001
t1 = period * dt_Integration
DPDtimes = f[:,0]*t1

plt.plot(DPDtimes,-f[:,1], "bo", markersize="0.5")
plt.axhline(y=0.0, color="r")

plt.xlabel('DPD times')
plt.ylabel('shear stress')

#plt.savefig('stress-strain.png',dpi=600, transparent=False)
plt.show()
