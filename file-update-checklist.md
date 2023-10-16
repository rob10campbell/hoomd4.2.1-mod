# Checklist for Updating Files

## Changes for adding flat and converging-diverging walls to the system
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] `methods/`
			* [ ] methods.py **wall amp freq, default amp=0 freq=0**
		* [ ] TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity, wall amp freq and bounceback**
		* [ ] TwoStepConstantVolume.h **wall amp and freq**


## adding shear rate for regular particles (does not include rigid bodies)
* [ ] `hoomd/`
	* [ ] BoxResizeUpdater.cc : **shear rate (SR)**
	* [ ] BoxResizeUpdater.h : **shear rate (SR), y-boundary velocity**
	* [ ] Communicator.cc : **shear rate (SR), y-boundary velocity**
	* [ ] Communicator.h : **shear rate (SR)**
	* [ ] ForceCompute.h : **shear rate (SR), virial_ind** 
	* [ ] HOOMDMPI.h : **shear rate (SR)** (add uint4)
	* [ ] Integrator.cc : **shear rate (SR)**
	* [ ] Integrator.h : **shear rate (SR)**
	* [ ] `md/`
		* [ ] EvaluatorPairMorse.h : **shear rate (SR)**
		* [ ] FIREEnergyMinimizer.cc : **shear rate (SR)**
		* [ ] FIREEnergyMinimizer.h : **shear rate (SR)**
		* [ ] integrate.py : **shear rate (SR), shear rate = 0 as default**
		* [ ] IntegrationMethodTwoStep.cc : **shear rate (SR)**
		* [ ] IntegrationMethodTwoStep.h : **shear rate (SR)**
		* [ ] IntegratorTwoStep.cc : **shear rate (SR)**
		* [ ] IntegratorTwoStep.h : **shear rate (SR)**
		* [ ] NeighborListTree.cc : **modifications for using our boundaries with MPI**
		* [ ] PotentialPairDPDThermo.h : **get box size/periodic, shear rate (SR)**
		* [ ] TwoStepBD.cc **shear rate (SR)**
		* [ ] TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity**
		* [ ] TwoStepLangevin.cc **shear rate (SR)**
	* [ ] `mpcd/`
		* [ ] Inegrator.cc : **shear rate (SR)**
		* [ ] Integrator.h : **shear rate (SR)**
	* [ ] `update/`
		* [ ] box_resize.py : **shear rate (SR)**
		
## tracking bond formation and breaking
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] CMakeLists.txt : **set new file (Lifetime.h)**
		* [ ] **[ADD NEW FILE]** Lifetime.h		
		* [ ] `pair`
			* [ ] pair.py : **call DPDMorse**
		* [ ] PotentialPairDPDThermo.h : **Lifetime, bond_calc, get particle diameter**

## Asakura-Oosawa Potential (might be incorrect calc?)
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] CMakeLists.txt : **set new file (EvaluatorPairDPDThermoDPDAO.h)**
		* [ ] **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDAO.h (the Asakura-Oosawa potential)
		* [ ] `pair`
			* [ ] \_\_init\_\_.py **call DPDAO**
			* [ ] pair.py : **call DPDAO**


## "core" mods (Contact force, Lubrication Force, track virial components)
* [ ] `hoomd/`
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
		* [x] module-md.cc : **add void/export for DPDMorse()**	
		* [x] `pair`
			* [x] \_\_init\_\_.py **call DPDMorse**
			* [x] pair.py : **define DPDMorse Python interface**
		* [x] PotentialPairDPDThermo.h : **virial_ind**
		* [ ] TwoStepBD.cc **contact force**
		* [ ] TwoStepLangevin.cc **contact force**
	* [x] ParticleData.cc : **virial_ind**
	* [x] ParticleData.h : **virial_ind**
