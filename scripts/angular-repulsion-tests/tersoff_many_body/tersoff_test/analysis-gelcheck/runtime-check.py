import math
import gsd.hoomd

filename = '../Gelation.gsd'

traj = gsd.hoomd.open(filename, 'rb')
nframes = len(traj)-1 # first frame is the 0th frame

# diffusion coefficients
kT = 0.1 # from simulation
R_C1 = 1 # radius of type 1 colloid (from simulation)
R_C2 = 2 # radius of type 2 colloid (from simulation)
eta0 = 0.3 # from simulation
d = 3 # dimension of the system (2D = 2, 3D = 3)

D_C1 = kT/(6*math.pi*eta0*R_C1) # diffusion coefficient = r^2 / 2d*tau
tau_C1 = (R_C1**2)/(2*d*D_C1)

D_C2 = kT/(6*math.pi*eta0*R_C2) # diffusion coefficient = r^2 / 2d*tau
tau_C2 = (R_C2**2)/(2*d*D_C2)

period = 10000
dt_Integration = 0.001
t1 = period * dt_Integration

sim_runtime_DPD = nframes*t1

sim_runtime_diff_C1 = nframes*t1/tau_C1
Brownian_times_100_C1 = 100*tau_C1
Btimes_timesteps_C1 = Brownian_times_100_C1 / dt_Integration

sim_runtime_diff_C2 = nframes*t1/tau_C2
Brownian_times_100_C2 = 100*tau_C2
Btimes_timesteps_C2 = Brownian_times_100_C2 / dt_Integration



print("----------")
print("1 Brownian time = 1 bare particle diffusion time (colloid 1) = "+str(round(tau_C1,2))+" DPD times")
print("100 Brownian times = 100 bare particle diffusion times (colloid 1) = "+str(round(Brownian_times_100_C1,2))+" DPD times")
print("In order to run for 100 Brownian times (colloid 1), sim needs to run for "+str(math.ceil(Btimes_timesteps_C1))+" timesteps\n")

print("1 Brownian time = 1 bare particle diffusion time (colloid 2) = "+str(round(tau_C2,2))+" DPD times")
print("100 Brownian times = 100 bare particle diffusion times (colloid 2) = "+str(round(Brownian_times_100_C2,2))+" DPD times")
print("In order to run for 100 Brownian times (colloid 2, sim needs to run for "+str(math.ceil(Btimes_timesteps_C2))+" timesteps\n")

print("sim ran for: \n  "+
     str(nframes*period)+" timesteps ; "+
     str(sim_runtime_DPD)+" DPD times \n  "+
     str(round(sim_runtime_diff_C1,2))+" diffiusion times (colloid 1) \n  "+
     str(round(sim_runtime_diff_C2,2))+" diffiusion times (colloid 2)")
print("----------")

if sim_runtime_DPD < Brownian_times_100_C1:
  remaining_time = math.ceil((Brownian_times_100_C1 - sim_runtime_DPD) * (1/dt_Integration))
  print("Needs more time. Run for an additional "+str(remaining_time)+" timesteps to reach 100 colloid 1 Brownian times.")
if sim_runtime_DPD < Brownian_times_100_C2:
  remaining_time = math.ceil((Brownian_times_100_C2 - sim_runtime_DPD) * (1/dt_Integration))
  print("Needs more time. Run for an additional "+str(remaining_time)+" timesteps to reach 100 colloid 2 Brownian times.")
elif (sim_runtime_DPD >= Brownian_times_100_C1) and (sim_runtime_DPD >= Brownian_times_100_C2):
  print("Simulation has run for at least 100 Brownian times. (colloid 1 AND colloid 2)")
print("----------")

