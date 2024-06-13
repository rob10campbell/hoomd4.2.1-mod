# From a GSD file, plot the shear stress (-pressure_xy) 
# of a simulation versus DPD time
# NOTE: requires "gsd-properties.txt" created with the
#	sim-analysis.py script

import numpy as np
import matplotlib.pyplot as plt

f = np.genfromtxt('gsd-properties.txt', skip_header=1)
plt.plot(f[:,0],-f[:,1], "bo", markersize="0.5")
plt.axhline(y=0.0, color="r")
plt.xlabel('Strains')
plt.ylabel('shear stress')
#plt.savefig('stres-strain.png',dpi=600, transparent=False)
plt.show()
