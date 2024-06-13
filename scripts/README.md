## Scripts

This folder contains example scripts for simulations in the PRO-CF Colloids Team.

1. [install-update](/scripts/install-update): bash scripts for installing and updating hoomd4.1-mod
2. [sim-templates](/scripts/sim-templates): template scripts for monodisperse colloidal simulations using [Brownian Dynamics](/scripts/sim-templates/Brownian-sims) or [Dissipative Particle Dynamics (DPD)](/scripts/sim-templates/DPD-sims). 
	* Options exist for initializing with a fixed simulation box size OR a fixed number of colloids. 
	* Initialization, equilibrium, and gelation scripts are provided for all methods. Shearing scripts are currently only available for DPD sims. 
	* Equilibrium sims include basic analysis scripts (extracting data from a GSD file, plotting temperature and shear)
	* Gelation sims include the most common advanced analyses that we use for particle data (more details on each of these is included in the analysis scripts):
		- extracting data from a GSD file (temperature and shear)
		- coordination number (Z), distribution and average value
		- mean squared displacement (MSD)
		- radial distribution function g(r)
		- qualitatively approximating the pair correlation function h(r1,r2)
		- the static structure factor S(q) (and associated g(r))
		- number density fluctuations (NDF), including options for various box sizes
		- pore size distribiution (PSD), AKA void size distribution
		- Voronoi volumes (distribution and statistical analysis)
		- largest connected component (LCC) and basic network science properties
3. [poly-templates](/scripts/poly-templates): template scripts for polydisperse colloidal simulations using [Brownian Dynamics](/scripts/poly-templates/BD-poly) or [Dissipative Particle Dynamics (DPD)](/scripts/poly-templates/DPD-poly) and a fixed-box. 
4. [shear-templates](/scripts/shear-templates): template scripts for shear simulations using Brownian Dynamics or Dissipative Particle Dynamics (DPD). 
	* Include analysis scripts for:
		- correcting DPD temperature and pressure values for colloids parunder shear
		- velocity profile
		- affinity of the motion (affine vs. non-affine)
		- fabric tensor
