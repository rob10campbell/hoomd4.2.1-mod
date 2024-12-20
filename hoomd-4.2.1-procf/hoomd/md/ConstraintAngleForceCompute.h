// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CONSTRAINTANGLEFORCECOMPUTE_H__
#define __CONSTRAINTANGLEFORCECOMPUTE_H__

#include "hoomd/BondedGroupData.h"
#include "hoomd/ForceCompute.h"
#include "NeighborList.h"

#include <memory>
#include <vector>
#include <set>
/*! \file ConstraintAngleForceCompute.h
    \brief Declares a class for computing harmonic angles
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace md
    {
struct angle_harmonic_params
    {
    Scalar k;
    Scalar t_0;

#ifndef __HIPCC__
    angle_harmonic_params() : k(0), t_0(0) { }

    angle_harmonic_params(pybind11::dict params)
        : k(params["k"].cast<Scalar>()), t_0(params["t0"].cast<Scalar>())
        {
        }

    pybind11::dict asDict()
        {
        pybind11::dict v;
        v["k"] = k;
        v["t0"] = t_0;
        return v;
        }
#endif
    }
#ifdef SINGLE_PRECISION
    __attribute__((aligned(8)));
#else
    __attribute__((aligned(16)));
#endif

//! Computes harmonic angle forces on each particle
/*! Constraint angle forces are computed on every particle in the simulation.

    The angles which forces are computed on are accessed from ParticleData::getAngleData
    \ingroup computes
*/
class PYBIND11_EXPORT ConstraintAngleForceCompute : public ForceCompute
    {
    public:
    //! Constructs the compute
    ConstraintAngleForceCompute(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<NeighborList> nlist);

    //! Destructor
    virtual ~ConstraintAngleForceCompute();

    //! Set the parameters
    virtual void setParams(unsigned int type, Scalar K, Scalar t_0);

    virtual void setParamsPython(std::string type, pybind11::dict params);

    /// Get the parameters for a type
    pybind11::dict getParams(std::string type);

#ifdef ENABLE_MPI
    //! Get ghost particle fields requested by this pair potential
    /*! \param timestep Current time step
     */
    virtual CommFlags getRequestedCommFlags(uint64_t timestep)
        {
        CommFlags flags = CommFlags(0);
        flags[comm_flag::tag] = 1;
        //flags[comm_flag::diameter] = 1;
        flags |= ForceCompute::getRequestedCommFlags(timestep);
        return flags;
        }
#endif

    protected:
    std::shared_ptr<NeighborList> m_nlist;
    Scalar* m_K;   //!< K parameter for multiple angle tyes
    Scalar* m_t_0; //!< r_0 parameter for multiple angle types

    std::shared_ptr<AngleData> m_angle_data; //!< Angle data to use in computing angles
    std::set< std::vector<unsigned int> > m_angle_members;
    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);
    };

    } // end namespace md
    } // end namespace hoomd

#endif
