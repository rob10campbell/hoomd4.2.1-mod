import hoomd
import hoomd.md
import hoomd.custom
import numpy as np
import math
import gsd
import gsd.pygsd
import gsd.hoomd
import random
import sys
import os
################
#For steady shear, you just need to change to constant variant.
################
Amp = 0.1 #float(sys.argv[1])
t_ramp = int(1e5)
omega = 2*np.pi/t_ramp
#omega = float(sys.argv[2])
D0 = 12 #float(sys.argv[3])
alpha = 30 #float(sys.argv[4])
dt_Integration = 0.0001#float(sys.argv[3])

kT=0.1; eta=1.0
R_C1 = 1.0
vf = 0.1
N_C1 = 100000
V_C1=(4./3.) * math.pi * R_C1 ** 3
V_total = N_C1 * V_C1/vf
L_X = V_total**(1./3.)
L_Y = L_X; L_Z = L_X

gamma = 6.0*np.pi*eta*R_C1 # BD stock friction coefficient

vinf = Amp * omega * L_Y
Iter_cycle = int(2*np.pi/omega / dt_Integration)
nframe2= int(Iter_cycle/100)

device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=50)
sim.timestep=0
if os.path.exists('Gelation1.gsd'): #('../u%d_k%d/Gelation1.gsd'%(D0,alpha)):
    sim.create_state_from_gsd(filename='Gelation1.gsd') #filename='../u%d_k%d/Gelation1.gsd'%(D0,alpha))
else:
    sim.create_state_from_gsd(filename='Gelation1.gsd') #filename='../u%d_k%d/Gelation.gsd'%(D0,alpha))

all_ = hoomd.filter.Type(['A'])

nl = hoomd.md.nlist.Tree(buffer=0.05);
morse = hoomd.md.pair.Morse(default_r_cut=1.0, nlist=nl)
morse.params[('A', 'A')] = dict(D0=D0*kT, alpha=alpha, r0=2*R_C1, f_contact=100)
morse.r_cut[('A', 'A')] = 2.0*R_C1 + (7.0/alpha)

brownian=hoomd.md.methods.Brownian(filter=all_, kT=kT, default_gamma=gamma)
#brownian=hoomd.md.methods.Brownian(filter=all_, kT=kT, alpha=3.0*np.pi*eta)
integrator=hoomd.md.Integrator(dt=dt_Integration,vinf=hoomd.variant.Cosinusoid(vinf,sim.timestep,omega*dt_Integration),forces=[morse],methods=[brownian])
sim.operations.integrator = integrator

logger = hoomd.logging.Logger()
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
sim.operations.computes.append(thermodynamic_properties)
logger.add(thermodynamic_properties,quantities=['kinetic_temperature','pressure_tensor','potential_energy'])
logger.add(sim,quantities=['tps'])

gsd_writer2 = hoomd.write.GSD(filename='Shearing_%.1f_%.3f_u%d_k%d.gsd'%(omega,Amp,D0,alpha),trigger=nframe2,filter=all_,dynamic=['property','momentum','attribute'],mode='wb')
sim.operations.writers.append(gsd_writer2)
gsd_writer2.log = logger

box_resize=hoomd.update.BoxShear(trigger=1,erate=hoomd.variant.Cosinusoid(omega*Amp,sim.timestep,omega*dt_Integration),deltaT=dt_Integration,flip=False)
sim.operations += box_resize
#sim.run(Iter_cycle*5)
n_cycles = 2
for i in range(n_cycles):
  sim.run(Iter_cycle)
  gsd_writer.flush()

