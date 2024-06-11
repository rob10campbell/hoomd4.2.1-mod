## Extract data from a GSD file to analyze the results 
## of a DPD colloid simulation
## NOTE: this code assumes 1 solvent type (typeid=0) 
##        and 1 colloid type (typeid=1)
##
## Extracts temperature and pressure (-stress) data from a GSD file 
## (Rob Campbell)


######### MODULE LIBRARY
import numpy as np
import gsd.hoomd
import math

######### INPUT PARAMETERS
filepath = '../Equilibrium-poly-DPD.gsd'

######### DEFINE SIM CHECKS
## extract thermodynamic properties 
## (temperature and compoentns of pressure (AKA -stress))
def extract_properties_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframe = len(traj)
	
  # create a file to save the thermodynamic properties
  f=open('gsd-properties.txt', 'w')
  f.write('simframe Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps\n')

  # for each frame
  for i in range(0, nframe):
    # extract the total "force virial contribution" of the pressure tensor
    Pr=float(traj[i].log['md/compute/ThermodynamicQuantities/pressure_tensor'][1])
    # extract the decomposed virial compoenents of the pressure tensor
    Vr=np.asarray(traj[i].log['md/compute/ThermodynamicQuantities/virial_ind_tensor'],dtype=float)
    # extract the potential energy and scale to kilo-units (1/1000)
    Pe=float(traj[i].log['md/compute/ThermodynamicQuantities/potential_energy'][0])*0.001
    # extract the kinetic temperature (kT)
    KT=float(traj[i].log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0])
    # extract the transactions per secont (tps) for assessing speed/efficiency of the simulation
    tps=int(traj[i].log['Simulation/tps'][0])

    # write these values to the output file:
    # raw values
    f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(i, Pr, 
      Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))
		
    # rounded values
    #f.write('%f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %d %0.2f %d\n'%(i, Pr, 
    #  Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))


######### RUN CHECKS ON A SIMULATION

if __name__ == '__main__':
  extract_properties_py(filepath)
