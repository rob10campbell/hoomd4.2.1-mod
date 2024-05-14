# Changelog with full list of changed files


If the box next to a filename is NOT checked, that means the modifications have not been added

File lists are formatted as: `folder/`; file

* [Core Mods](/changelog.md#core-modifications)
* [Shear Rate DPD](/changelog.md#shear-rate-DPD)
* [Shear Rate BD](/changelog.md#shear-rate-BD)
* [Bond tracking](/changelog.md#bond-tracking)
* [Walls](/changelog.md#walls)
* [Asakura-Oosawa Potential](/changelog.md#asakura-oosawa-potential)

## Core Modifications
Contact Force, Lubrication Force, track virial components (Nabi and Deepak)
- **DPDMorse**: Add a new method for evaluating the pair-forces between two particles called "DPDMorse." This method calculates the correct combination of forces for each pair of particles as described in the [background reading on DPD for Colloids](/background-reading/2-DPD-for-Colloids-18pg.pdf); i.e., the standard DPD forces plus the Morse potential, hydrodynamics (the squeezing force AKA lubrication force), and a contact force for resolving semi-hard colloid-colloid particle overlaps. (filename: `EvaluatorPairDPDThermoDPDMorse.h`)
- **contact force**: Add a contact force to Brownian simulations for resolving semi-hard colloid-colloid particle overlaps
- **virial_ind**: Add the ability to track the "independent virials" (AKA the virial_ind) for each particle pair. This is the contribution of each of the individual forces in the virial (conservative, dissipative, random, Morse, lubrication, and contact) in addition to the total virial component of the pressure tensor

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
			* [x] pair.py : **define DPDMorse Python interface, add Morse contact force param**
		* [x] PotentialPairDPDThermo.h : **virial_ind**
	* [x] ParticleData.cc : **virial_ind**
	* [x] ParticleData.h : **virial_ind**

## Shear Rate DPD
Adding shear rate for regular particles (does nor include rigid bodies) (Deepak)
- **shear rate (SR)**: Add a new class for controlling the shear rate (SR) when shearing a system for DPDMorse and Morse Brownian sims (BD and Langevin integration); however this does not apply to rigid bodies
- **Tree modifications**: Modify the Tree Neighboring List method to remove auto-updates without SR
- **y-boundary velocity**: Modify the way that a particle's velocity is updated when it crosses a y-boundary so that it takes into account the effect of the applied shear rate (as described in the [background reading on shearing](/background-reading/4-Shearing-4pg.pdf)

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
			* [x] test_communication.cc **shear rate (SR)**
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

## Shear Rate BD
Adding shear rate for regular particles
* [ ] `hoomd/`
	* [ ] BoxShearUpdater.cc
	* [ ] BoxShearUpdater.h
	* [ ] CmakeList.txt
 	* [ ] Variant.cc
	* [ ] Variant.h
	* [ ] `update/`
 		  * [ ] init.py (add BoxShear?)
	 	  * [ ] box_shear.py
	* [ ] variant.py

## Bond Tracking
Track bond formation and breaking (Nabi, Deepak, and Rob)
- **Lifetime**: Add the ability to track the formation and breaking of any bond between two particles, and calculate and record the bond lifetime (filename: `Lifetime.h`)
- **bond_calc**: Create a flag to turn the bond-lifetime calculation On or Off (True/False) *NOTE: lifetime calculation is slow and was not able to be sped up with MPI*

* [x] `hoomd/`
	* [x] `md/`
		* [x] CMakeLists.txt : **set new file (Lifetime.h)**
		* [x] **[ADD NEW FILE]** Lifetime.h		
		* [x] `pair`
			* [x] pair.py : **add bond_calc flag to Pair and DPDMorse**
		* [x] PotentialPairDPDThermo.h : **Lifetime, bond_calc, get particle diameter (for bond calc)**

## Walls
Wall options: flat or converging diverging (Josh)
- **walls**: Add walls in the Y direction that are either flat or sinusoidal in shape, including correct bounce-back conditions for solvents interacting with the wall. Add the ability for opposite sinusoidal walls to be aligned (in phase) or offset by a phase shift and for the "sides" to have flat walls.
* [ ] `hoomd/`
	* [ ] `md/`
		* [ ] `methods/`
			* [ ] methods.py **wall amp freq, default amp=0 freq=0**
		* [ ] TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity, wall amp freq and bounceback**
		* [ ] TwoStepConstantVolume.h **wall amp and freq**

		
## Asakura-Oosawa Potential 
Adding AO Potential (might be incorrect calc?) (Rob)
- **DPDAO**: Add the Asakura-Oosawa potential as an alternative to Morse potential. This is an laternate method for evaluating the pair-forces between two particles that calculates the correct combination of forces for each pair of particles as described in the [background reading on DPD for Colloids](/background-reading/2-DPD-for-Colloids-18pg.pdf), but replaces the Morse Potential with the Asakura-Oosawa (AO) potential; i.e., the standard DPD forces plus the AO potential, hydrodynamics (the squeezing force AKA lubrication force), and a contact force for resolving semi-hard colloid-colloid particle overlaps.
 * [ ] `hoomd/`
	* [ ] `md/`
		* [ ] CMakeLists.txt : **set new file (EvaluatorPairDPDThermoDPDAO.h)**
		* [ ] **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDAO.h (the Asakura-Oosawa potential)
		* [ ] `pair`
			* [ ] \_\_init\_\_.py **call DPDAO**
			* [ ] pair.py : **call DPDAO**
