# hoomd4.2.1-mod

This repository contains the modified version of HOOMD-blue v4.2 used for colloid simulations in the PRO-CF group. It also includes: 
* [Installation instructions](/README.md#installation)
* [Background Reading](/background-reading) about how these simulations work
* [Example scripts](/scripts) for installing and running simulations on an HPC research cluster
* [Citation information](/citation-guide.md) for papers published using this simulation platform
* A [Changelog](/README.md#changelog) summarizing what was changed
* [A full list of the files that were changed](/README.md#Full-list-of-changed-files)

Additional branches are available, tho they may be incomplete:
- branch "polydispersity": modifications for non-uniform particle sizes
- branch "hoomd4.2_w_wall": modifications for flat and sinusoidal walls
- branch "no_shear": clean version without shear or bond-tracking

[Last Updated: December 2023]

Contact: Rob Campbell (campbell.r@northeastern.edu)

-----------------
For any questions, or to help with modifications, contact Rob.

To-Do:
- [x] copy in core mods 
- [x] test compile
- [x] resolve segmentation fault in DPDMorse or Integrator
- [x] update sim template scripts
- [x] test initialization, equilibrium, and gelation with DPDMorse
- [x] test initialization, equilibrium, and gelation with BD/Langevin
- [x] copy in Lifetime.h
- [x] test compile
- [x] resolve Lifetime compilation error re: "getMPICommunicator()"
- [x] test gelation sim with bond lifetime tracking
- [ ] add bond tracking to BD/Langevin sims?
- [x] copy in shear mods
- [x] test compile
- [x] test shear sim with DPDMorse
- [ ] shear does not work in BD/Langevin (thermostat issues)
- [ ] copy in wall mods
- [ ] test compile
- [ ] test DPD initialization, equilibrium, and gelation w/ and w/out walls
- [ ] test BD/Langevin sims with and without walls
- [ ] remove AO from README or add AO mods
- [ ] update sim analysis scripts
- [ ] test sim analysis
- [x] update colloids-setup repo to use hoomd4.2 (!)
-----------------

## Installation

To install hoomd4.2.1-mod on Discovery:

Login and move to the location you want to put hoomd4.2-mod
```bash
ssh your-username@login.discovery.neu.edu
```
```bash
cd /work/props/your-name/software/
```
Make a new directory for hoomd4.2.1-mod
```bash
mkdir hoomd4.2.1-mod
```
Clone the repository
```bash
git clone git@github.com:procf/hoomd4.2.1-mod.git
```
Run the install-hoomd-mod-4.2.1 script with sbatch
```bash
cd /scripts/install-update/ && sbatch install-hoomd-mod-4.2.1
```
<br>

## Changelog

Core Modifications: Contact Force, Lubrication Force, track virial components (Nabi and Deepak)
- **DPDMorse**: Add a new method for evaluating the pair-forces between two particles called "DPDMorse." This method calculates the correct combination of forces for each pair of particles as described in the [background reading on DPD for Colloids](/background-reading/2-DPD-for-Colloids-18pg.pdf); i.e., the standard DPD forces plus the Morse potential, hydrodynamics (the squeezing force AKA lubrication force), and a contact force for resolving semi-hard colloid-colloid particle overlaps. (filename: `EvaluatorPairDPDThermoDPDMorse.h`)
- **virial_ind**: Add the ability to track the "independent virials" (AKA the virial_ind) for each particle pair. This is the contribution of each of the individual forces in the virial (conservative, dissipative, random, Morse, lubrication, and contact) in addition to the total virial component of the pressure tensor

Track Bond Formation and Breaking (Nabi, Deepak, and Rob)
- **Lifetime**: Add the ability to track the formation and breaking of any bond between two particles, and calculate and record the bond lifetime (filename: `Lifetime.h`)
- **bond_calc**: Create a flag to turn the bond-lifetime calculation On or Off (True/False) *NOTE: lifetime calculation is slow and was not able to be sped up with MPI*

Shear Rate Modification (Deepak)
- **shear rate (SR)**: Add a new class for controlling the shear rate (SR) when shearing a system for DPDMorse and Morse Brownian sims (BD and Langevin integration); however this does not apply to rigid bodies
- **Tree**: Modify the Tree Neighboring List method to remove auto-updates without SR
- **y-boundary velocity**: Modify the way that a particle's velocity is updated when it crosses a y-boundary so that it takes into account the effect of the applied shear rate (as described in the [background reading on shearing](/background-reading/4-Shearing-4pg.pdf)

Wall options: flat or converging diverging (Josh)
- **walls**: Add walls in the Y direction that are either flat or sinusoidal in shape, including correct bounce-back conditions for solvents interacting with the wall. Add the ability for opposite sinusoidal walls to be aligned (in phase) or offset by a phase shift and for the "sides" to have flat walls.

Asakura-Oosawa Potential (Rob)
- **DPDAO**: Add the Asakura-Oosawa potential as an alternative to Morse potential. This is an laternate method for evaluating the pair-forces between two particles that calculates the correct combination of forces for each pair of particles as described in the [background reading on DPD for Colloids](/background-reading/2-DPD-for-Colloids-18pg.pdf), but replaces the Morse Potential with the Asakura-Oosawa (AO) potential; i.e., the standard DPD forces plus the AO potential, hydrodynamics (the squeezing force AKA lubrication force), and a contact force for resolving semi-hard colloid-colloid particle overlaps.
 

For more details on which files were changes to accomodate these changes, see the [full list of changed files](#full-list-of-changed-files) at the end of this document.
<br>
<br>
## Full list of changed files

(formatted as: `folder/`; file)

* `hoomd/`
	* BoxResizeUpdater.cc : **shear rate (SR)**
	* BoxResizeUpdater.h : **shear rate (SR), y-boundary velocity**
	* Communicator.cc : **shear rate (SR), virial_ind, y-boundary velocity**
	* Communicator.h : **virial_ind, shear rate (SR)**
	* ForceCompute.cc : **virial_ind**
	* ForceCompute.h : **shear rate (SR), virial_ind** 
	* HOOMDMPI.h : **shear rate (SR)** (add uint4)
	* Integrator.cc : **shear rate (SR), virial_ind**
	* Integrator.h : **shear rate (SR)**
	* `md/`
		* CMakeLists.txt : **set new files (EvaluatorPairDPDThermoDPDMorse.h, EvaluatorPairDPDThermoDPDAO.h, and Lifetime.h)**
		* compute.py : **virial_ind**
		* ComputeThermo.cci : **virial_ind**
		* ComputeThermo.h : **virial_ind**
		* ComputeThermoTypes.h : **virial_ind**
		* **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDAO.h (the Asakura-Oosawa potential)
		* EvaluatorPairDPDThermoDPD.h : **virial_ind (redefine forces to include all virial_ind terms)**
		* **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDMorse.h
		* EvaluatorPairDPDThermoLJ.h : **virial_ind**
		* EvaluatorPairMorse.h : **shear rate (SR)**
		* FIREEnergyMinimizer.cc : **shear rate (SR)**
		* FIREEnergyMinimizer.h : **shear rate (SR)**
		* integrate.py : **shear rate (SR), shear rate = 0 as default**
		* IntegrationMethodTwoStep.cc : **shear rate (SR)**
		* IntegrationMethodTwoStep.h : **shear rate (SR)**
		* IntegratorTwoStep.cc : **shear rate (SR)**
		* IntegratorTwoStep.h : **shear rate (SR)**
		* **[ADD NEW FILE]** Lifetime.h
		* module-md.cc : **add void/export for DPDMorse() and DPDAO()**
		* NeighborListTree.cc : **modifications for using our boundaries with MPI**
		* `pair`
			* \_\_init\_\_.py **call DPDMorse and DPDAO**
			* pair.py : **call DPDMorse and DPDAO**
		* PotentialPairDPDThermo.h : **Lifetime, bond_calc, get particle diameter, virial_ind, get box size/periodic, shear rate (SR)**
  		* `test/`
    			* test_fire_energy_minimizer.cc **shear rate (SR)** 
		* TwoStepBD.cc **shear rate (SR)**
		* TwoStepConstantVolume.cc : **shear rate (SR), y-boundary velocity**	
		* TwoStepLangevin.cc **shear rate (SR)**
	* `mpcd/`
		* Inegrator.cc : **shear rate (SR)**
		* Integrator.h : **shear rate (SR)**
	* ParticleData.cc : **virial_ind**
	* ParticleData.h : **virial_ind**
  	* `test/`
        	* test_system.cc **shear rate (SR)** 
	* `update/`
		* box_resize.py : **shear rate (SR)**
