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
	* [ ] Communicator.cc : **virial_ind**
	* [ ] Communicator.h : **virial_ind**
	* [ ] ForceCompute.cc : **virial_ind**
	* [ ] ForceCompute.h : **virial_ind** 
	* [ ] HOOMDMPI.h : ? **shear rate (SR)** (add uint4)
	* [ ] Integrator.cc : **virial_ind**	
	* [ ] `md/`
		* [ ] CMakeLists.txt : **set new file (EvaluatorPairDPDThermoDPDMorse.h)**
		* [ ] compute.py : **virial_ind**
		* [ ] ComputeThermo.cc : **virial_ind**
		* [ ] ComputeThermo.h : **virial_ind**
		* [ ] ComputeThermoTypes.h : **virial_ind**
		* [ ] EvaluatorPairDPDThermoDPD.h : **virial_ind (redefine forces to include all virial_ind terms)**
		* [ ] **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDMorse.h
		* [ ] EvaluatorPairDPDThermoLJ.h : **virial_ind**
		* [ ] module-md.cc : **add void/export for DPDMorse()**	
		* [ ] `pair`
			* [ ] \_\_init\_\_.py **call DPDMorse**
			* [ ] pair.py : **call DPDMorse**
		* [ ] PotentialPairDPDThermo.h : **get particle diameter, virial_ind**
	* [ ] ParticleData.cc : **virial_ind**
	* [ ] ParticleData.h : **virial_ind**
