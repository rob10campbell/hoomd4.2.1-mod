// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef _COMPUTE_THERMO_TYPES_H_
#define _COMPUTE_THERMO_TYPES_H_

#include "hoomd/HOOMDMath.h"
/*! \file ComputeThermoTypes.h
    \brief Data structures common to both CPU and GPU implementations of ComputeThermo
    */

namespace hoomd
    {
namespace md
    {
//! Enum for indexing the GPUArray of computed values
struct thermo_index
    {
    //! The enum
    enum Enum
        {
        translational_kinetic_energy = 0, //!< Index for the kinetic energy in the GPUArray
        rotational_kinetic_energy,        //!< Rotational kinetic energy
        potential_energy,                 //!< Index for the potential energy in the GPUArray
        pressure,                         //!< Total pressure
        pressure_xx,   //!< Index for the xx component of the pressure tensor in the GPUArray
        pressure_xy,   //!< Index for the xy component of the pressure tensor in the GPUArray
        pressure_xz,   //!< Index for the xz component of the pressure tensor in the GPUArray
        pressure_yy,   //!< Index for the yy component of the pressure tensor in the GPUArray
        pressure_yz,   //!< Index for the yz component of the pressure tensor in the GPUArray
        pressure_zz,   //!< Index for the zz component of the pressure tensor in the GPUArray
        //~ add virial_ind [RHEOINF] 
	virial_ind_xx, //~!< Index for the conservative component of the virial_ind tensor in the GPUArray
        virial_ind_xy, //~!< Index for the dissipative component of the virial_ind tensor in the GPUArray
        virial_ind_xz, //~!< Index for the random component of the virial_ind tensor in the GPUArray
        virial_ind_yy, //~!< Index for the squeezing_hydro component of the virial_ind tensor in the GPUArray
        virial_ind_yz, //~!< Index for the contact component of the virial_ind tensor in the GPUArray
	//~

        num_quantities // final element to count number of quantities
        };
    };

//! structure for storing the components of the pressure tensor
struct PressureTensor
    {
    //! The six components of the upper triangular pressure tensor
    Scalar xx; //!< xx component
    Scalar xy; //!< xy component
    Scalar xz; //!< xz component
    Scalar yy; //!< yy component
    Scalar yz; //!< yz component
    Scalar zz; //!< zz component
    };

//~ add virial_ind [RHEOINF]
//! structure for storing the components of the virial_ind tensor
struct VirialTensor
    {
    //! The five components of the virial_ind tensor
    Scalar xx; //~ conservative
    Scalar xy; //~ dissipative
    Scalar xz; //~ random
    Scalar yy; //~ squeezing_hydro
    Scalar yz; //~ contact
    };
//~

    } // end namespace md
    } // end namespace hoomd
#endif
