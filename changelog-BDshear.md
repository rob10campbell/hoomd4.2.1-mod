# Changelog with full list of files changed for BD-shear

test change

If the box next to a filename is NOT checked, that means the modifications have not been added

File lists are formatted as: `folder/`; file

* [Shear Rate BD](/changelog.md#shear-rate-BD)

## virial_xyi_ind BD
- **charge** make sure charge is available for body force
* [x] `hoomd/`
	* [x] `md/`
		* [x] PotentialPairDPDThermo.h : **charge**


## Shear Rate BD
New way to add shear rate for regular particles, that works for both DPD and BD (does not include rigid bodies) (Deepak)
- **flow velocity** add the flow velocity for calculating shear rate
- **y-boundary velocity** update velocity when wrapping particles across the Y boundary (gradient direction)
- **shear rate** add the shear rate variable
- **iostream** import package for Variant class
- **Oscillatory** create variant class for t_ramp oscillatory shear (legacy)
- **Sinusoidal** create variant class for sinusoidal oscillation (legacy)
- **Cosinusoidal** create variant class for cosinusoidal oscillation (use this one)
- **BoxShear** and **box_shear.py** new classes for shearing
- **virialxyi_ind** track xyi component of virial separately in Brownian shear sims


* [x] `hoomd/`
	* [x] BoxResizeUpdater.h : **flow velocity (vinf, cur_vel)**
	* [x] BoxResizeUpdater.cc : **flow velocity (vinf, cur_vel), y-boundary velocity**
	* [x] **[ADD NEW FILE]** BoxShearUpdater.h
	* [x] **[ADD NEW FILE]** BoxShearUpdater.cc
	* [x] Communicator.cc : **shear rate (SR), y-boundary velocity**
	* [x] Communicator.h : **shear rate (SR)**
	* [x] ForceCompute.h : **shear rate (SR)**
	* [x] HOOMDMPI.h : **shear rate (SR)** (add uint4) 
	* [x] CmakeList.txt : **add BoxShearUpdater.h, BoxShearUpdater.cc**
	* [x] Integrator.h : **flow velocity (vinf), shear rate (SR)**
	* [x] Integrator.cc : **flow velocity (vinf), shear rate (SR), Box_Dim**
	* [x] module.cc : **BoxShear**
	* [x] `md/`
		* [x] FIREEnergyMinimizer.cc : **flow velocity (vinf)**
		* [x] FIREEnergyMinimizer.h : **flow velocity (vinf)**
		* [x] IntegrationMethodTwoStep.h : **shear rate (SR)**
		* [x] IntegrationMethodTwoStep.cc : **shear rate (SR)**
		* [x] IntegratorTwoStep.cc : **flow velocity (vinf)**
		* [x] IntegratorTwoStep.h : **flow velocity (vinf), shear rate (SR)**
		* [x] integrate.py : **flow velocity (vinf) and default to 0**
                * [x] PotentialPair.h : **virialxyi_ind**
		* [x] PotentialPairDPDThermo.h **shear rate (SR), virialxyi_ind**
		* [x] `test/`
			* [x] test_communication.cc : **dummy flow velocity (vinf)** 
			* [x] test_fire_energy_minimizer.cc : **dummy flow velocity (vinf)**
		* [x] TwoStepBD.cc : **flow velocity (vinf, cur_vel), y-boundary velocity**
		* [x] TwoStepConstantVolume.cc : **box dims, shear rate (SR)**
		* [x] TwoStepLangevin.cc : **flow velocity (vinf), shear rate (SR)**
	* [x] `mpcd/`
		* [x] Integrator.h : **flow velocity (vinf)**
		* [x] Integrator.cc : **flow velocity (vinf)**
	* [x] `test/` 
		* [x] test_system.cc : **dummy flow velocity (vinf)**
	* [x] `update/`
		* [x] CmakeList.txt : **add box_shear.py**
		* [x] \_\_init\_\_.py : **add BoxShear**
		* [x] **[ADD NEW FILE]** box_shear.py
		* [x] box_resize.py : **flow velocity (vinf)**
	* [x] Variant.h : **iostream, Oscillatory, Sinusoidal, Cosinusoidal**
	* [x] Variant.cc : **Oscillatory, Sinusoidal, Cosinusoidal**
	* [x] variant.py : **Oscillatory, Sinusoidal, Cosinusoidal**
