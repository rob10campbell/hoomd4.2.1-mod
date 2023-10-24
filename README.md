# hoomd4.2.1-mod

This repository contains the core modified version of HOOMD-blue v4.2 used for colloid simulations in the PRO-CF group. It also includes: 
* [Installation instructions](/README.md#installation)
* [Background Reading](/background-reading) about how these simulations work
* [Example scripts](/scripts) for installing and running simulations on an HPC research cluster
* [Citation information](/citation-guide.md) for papers published using this simulation platform
* A [Changelog](/README.md#changelog) summarizing what was changed
* [A full list of the files that were changed](/README.md#Full-list-of-changed-files)

[Last Updated: September 2023]

Contact: Rob Campbell (campbell.r@northeastern.edu)

To-Do:
- [ ] update sim template scripts
- [ ] test initialization, equilibrium, and gelation with BD/Langevin


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
**contact force**: add contact force for Brownian Dynamics

For more details on which files were changes to accomodate these changes, see the [full list of changed files](#full-list-of-changed-files) at the end of this document.
<br>
<br>
## Full list of changed files

(formatted as: `folder/`; file)

* `hoomd/`
	* Communicator.cc : **virial_ind**
	* Communicator.h : **virial_ind**
	* ForceCompute.cc : **virial_ind**
	* ForceCompute.h : **virial_ind** 
	* Integrator.cc : **virial_ind**
	* `md/`
		* CMakeLists.txt : **set new files (EvaluatorPairDPDThermoDPDMorse.h)**
		* compute.py : **virial_ind**
		* ComputeThermo.cc : **virial_ind**
		* ComputeThermo.h : **virial_ind**
		* ComputeThermoTypes.h : **virial_ind**
		* EvaluatorPairDPDThermoDPD.h : **virial_ind (redefine forces to include all virial_ind terms)**
		* **[ADD NEW FILE]** EvaluatorPairDPDThermoDPDMorse.h
		* EvaluatorPairDPDThermoLJ.h : **virial_ind**
		* EvaluatorPairMorse.h : **contact force**
		* module-md.cc : **add void/export for DPDMorse()**
		* `pair`
			* \_\_init\_\_.py **call DPDMorse**
			* pair.py : **call DPDMorse**
		* PotentialPairDPDThermo.h : **virial_ind**
	* ParticleData.cc : **virial_ind**
	* ParticleData.h : **virial_ind**
