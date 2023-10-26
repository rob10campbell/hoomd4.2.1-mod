# Checklist for Updating Files

## Changes for adding flat and converging-diverging walls to the system
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] `methods/`
			* [ ] methods.py **wall amp freq, default amp=0 freq=0**
		* [ ] TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity, wall amp freq and bounceback**
		* [ ] TwoStepConstantVolume.h **wall amp and freq**

## adding shear rate for regular particles (does not include rigid bodies)
* [x] `hoomd/`
	* [x] BoxResizeUpdater.cc : **shear rate (SR)**
	* [x] BoxResizeUpdater.h : **shear rate (SR), y-boundary velocity**
	* [x] Communicator.cc : **shear rate (SR), y-boundary velocity**
	* [x] Communicator.h : **shear rate (SR)**
	* [x] ForceCompute.h : **shear rate (SR), virial_ind** 
	* [x] HOOMDMPI.h : **shear rate (SR)** (add uint4)
	* [x] Integrator.cc : **shear rate (SR)**
	* [x] Integrator.h : **shear rate (SR)**
	* [x] `md/`
		* [x] FIREEnergyMinimizer.cc : **shear rate (SR)**
		* [x] FIREEnergyMinimizer.h : **shear rate (SR)**
		* [x] integrate.py : **shear rate (SR), shear rate = 0 as default**
		* [x] IntegrationMethodTwoStep.cc : **shear rate (SR)**
		* [x] IntegrationMethodTwoStep.h : **shear rate (SR)**
		* [x] IntegratorTwoStep.cc : **shear rate (SR)**
		* [x] IntegratorTwoStep.h : **shear rate (SR)**
		* [x] NeighborListTree.cc : **modifications for using our boundaries with MPI**
		* [x] PotentialPairDPDThermo.h : **get box size/periodic, shear rate (SR)**
  		* [x] `test/`
    			* [x] test_fire_energy_minimizer.cc **shear rate (SR)** 
		* [x] TwoStepBD.cc **shear rate (SR)**
		* [x] TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity**
		* [x] TwoStepLangevin.cc **shear rate (SR)**
	* [x] `mpcd/`
		* [x] Integrator.cc : **shear rate (SR)**
		* [x] Integrator.h : **shear rate (SR)**
  	* [x] `test/`
        	* [x] test_system.cc **shear rate (SR)** 
	* [x] `update/`
		* [x] box_resize.py : **shear rate (SR)**
		
## tracking bond formation and breaking
* [x] `hoomd/`
	* [x] `md/`
		* [x] CMakeLists.txt : **set new file (Lifetime.h)**
		* [x] **[ADD NEW FILE]** Lifetime.h		
		* [x] `pair`
			* [x] pair.py : **add bond_calc flag**
		* [x] PotentialPairDPDThermo.h : **Lifetime, bond_calc, get particle diameter (for bond calc)**

## Asakura-Oosawa Potential (might be incorrect calc?)
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] CMakeLists.txt : **set new file (EvaluatorPairDPDThermoDPDAO.h)**
		* [ ] **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDAO.h (the Asakura-Oosawa potential)
		* [ ] `pair`
			* [ ] \_\_init\_\_.py **call DPDAO**
			* [ ] pair.py : **call DPDAO**


## "core" mods (Contact force, Lubrication Force, track virial components)
* [x] `hoomd/`
	* [x] Communicator.cc : **virial_ind**
	* [x] Communicator.h : **virial_ind**
	* [x] ForceCompute.cc : **virial_ind**
	* [x] ForceCompute.h : **virial_ind** 
	* [x] Integrator.cc : **virial_ind**	
	* [x] `md/`
		* [x] CMakeLists.txt : **set new file (EvaluatorPairDPDThermoDPDMorse.h)**
		* [x] compute.py : **virial_ind**
		* [x] ComputeThermo.cc : **virial_ind**
		* [x] ComputeThermo.h : **virial_ind**
		* [x] ComputeThermoTypes.h : **virial_ind**
		* [x] EvaluatorPairDPDThermoDPD.h : **virial_ind (redefine forces to include all virial_ind terms)**
		* [x] **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDMorse.h
		* [x] EvaluatorPairDPDThermoLJ.h : **virial_ind**
		* [x] EvaluatorPairMorse.h : **contact force**
		* [x] module-md.cc : **add void/export for DPDMorse()**	
		* [x] `pair`
			* [x] \_\_init\_\_.py **call DPDMorse**
			* [x] pair.py : **define DPDMorse Python interface**
		* [x] PotentialPairDPDThermo.h : **virial_ind**
	* [x] ParticleData.cc : **virial_ind**
	* [x] ParticleData.h : **virial_ind**
