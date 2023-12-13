# hoomd4.2.1-mod

This repository contains the modified version of HOOMD-blue v4.2 used for colloid simulations in the PRO-CF group. It also includes: 
* [Installation instructions](/README.md#installation)
* [Background Reading](/background-reading) about how these simulations work
* [Example scripts](/scripts) for installing and running simulations on an HPC research cluster
* [Citation information](/citation-guide.md) for papers published using this simulation platform
* A [Changelog](/changelog.md) summarizing what was changed, and a full list of the files that were changed

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
