# hoomd4.2.1-mod

This repository contains the modified version of HOOMD-blue v4.2 used for colloid simulations in the PRO-CF group. It also includes: 
* [Installation instructions](/README.md#installation)
* [Background Reading](/background-reading) about how these simulations work
* [Example scripts](/scripts) for installing and running simulations on an HPC research cluster
* [Citation information](/citation-guide.md) for papers published using this simulation platform
* A [Changelog](/changelog.md) summarizing what was changed, and a full list of the files that were changed

Additional branches are available, tho they may be incomplete:
- branch "hoomd4.2_w_wall": modifications for flat and sinusoidal walls
- branch "no_shear": clean version of original DPDMorse and virial_ind tracking without shear, bond-tracking, or other mods.

[Last Updated: June 2024]

Contact: Rob Campbell (campbell.r@northeastern.edu)

-----------------
For any questions, or to help with modifications, contact Rob.

To-Do:
- [ ] fix bond trackign to work correctly @Rob
	- [ ] and add bond tracking to BD/Langevin sims?
- [ ] copy in wall mods @Josh @Rob
	- [ ] test compile wall mods
	- [ ] test DPD initialization, equilibrium, and gelation w/ walls
	- [ ] test BD/Langevin sims w/ walls
- [ ] add bond contraints and associated bond tracking @Paniz
- [ ] remove AO from README or add AO mods @Rob
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
