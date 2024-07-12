// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_OPP_H__
#define __PAIR_EVALUATOR_OPP_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h" // add vectors for optional position [RHEOINF]

/*! \file EvaluatorPairOPP.h
    \brief Defines the pair evaluator class for OPP potential
*/

// need to declare these class methods with __device__ qualifiers when building
// in nvcc DEVICE is __host__ __device__ when included in nvcc and blank when
// included into the host compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hoomd
    {
namespace md
    {
//! Class for evaluating the oscillating pair potential
/*! <b>General Overview</b>

    <b>OPP specifics</b>

    EvaluatorPairOPP evaluates the function:
    \f{equation*}
    V_{\mathrm{OPP}}(r)  = - C_1 r^{\eta_1} + C_2 r^{\eta_2}
                             \cos(k ( r - b) - \phi)
    \f}

*/
class EvaluatorPairOPP
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar C1;
        Scalar C2;
        Scalar eta1;
        Scalar eta2;
        Scalar k;
        Scalar phi;

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : C1(0), C2(0), eta1(0), eta2(0), k(0), phi(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            C1 = v["C1"].cast<Scalar>();
            C2 = v["C2"].cast<Scalar>();
            eta1 = v["eta1"].cast<Scalar>();
            eta2 = v["eta2"].cast<Scalar>();
            k = v["k"].cast<Scalar>();
            phi = v["phi"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["C1"] = C1;
            v["C2"] = C2;
            v["eta1"] = eta1;
            v["eta2"] = eta2;
            v["k"] = k;
            v["phi"] = phi;
            return v;
            }
#endif
        };

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _radcontact the sum of the interacting particle radii [RHEOINF]
        \param _pair_typeids the typeIDs of the interacting particles [RHEOINF]
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairOPP(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~ add radcontact, pair_typeIDs [RHEOINF]
        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), params(_params) //~ add radcontact [RHEOINF]
        {
        typei = _pair_typeids[0]; //~ add typei [RHEOINF]
        typej = _pair_typeids[1]; //~ add typej [RHEOINF]
        }

    //~ add tags [RHEOINF]
    HOSTDEVICE static bool needsTags()
        {
        return false;
        }
    HOSTDEVICE void setTags(unsigned int tagi, unsigned int tagj) { }
    //~

    //!~ add diameter [RHEOINF]
    DEVICE static bool needsDiameter()
        {
        return false;
        }
    //! Accept the optional diameter values
    /*! \param di Diameter of particle i
        \param dj Diameter of particle j
    */
    DEVICE void setDiameter(Scalar di, Scalar dj) { }
    //~

    //! OPP doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }

    //! Accept the optional charge values.
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //~ add timestep [RHEOINF]
    HOSTDEVICE static bool needsTimestep()
        {
        return false;
        }
    HOSTDEVICE void setTimestep(uint64_t timestep)  { }
    //~

    //~ add i and j positions [RHEOINF] 
    DEVICE static bool needsIJPos()
        {
        return false;
        }
    //! Accept the optional position values
    /*! \param pi position of particle i
        \param pj position of particle j
    */
    DEVICE void setIJPos(Scalar3 pi, Scalar3 pj) {}
    //~

    //!~ Whether the potential pair needs BoxDim info [RHEOINF]
    HOSTDEVICE static bool needsBox()
        {
        return false;
        }
    //! Accept the optional BoxDim structure
    /*! \param box the current box
    */
    HOSTDEVICE void setBox(const BoxDim box) { }
    //~

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force
     * divided by r.
     *  \param pair_eng Output parameter to write the computed pair energy
     *  \param energy_shift If true, the potential must be shifted so that
     *      V(r) is continuous at the cutoff

     *  \return True if they are evaluated or false if they are not because
     *  we are beyond the cutoff
     */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        if (rsq < rcutsq)
            {
            // Get quantities need for both energy and force calculation
            Scalar r(fast::sqrt(rsq));
            Scalar eval_sin, eval_cos;
            fast::sincos(params.k * r - params.phi, eval_sin, eval_cos);

            // Compute energy
            Scalar r_eta1_arg(params.C1 * fast::pow(r, -params.eta1));
            Scalar r_to_eta2(fast::pow(r, -params.eta2));
            Scalar r_eta2_arg(params.C2 * r_to_eta2 * eval_cos);
            pair_eng = r_eta1_arg + r_eta2_arg;

            // Compute force
            force_divr = ((r_eta1_arg * params.eta1 + r_eta2_arg * params.eta2) / r
                          + (params.C2 * params.k * r_to_eta2 * eval_sin))
                         / r;

            if (energy_shift)
                {
                Scalar r_cut(fast::sqrt(rcutsq));
                Scalar r_cut_eta1_arg(params.C1 * fast::pow(r_cut, -params.eta1));
                Scalar r_cut_eta2_arg(params.C2 * fast::pow(r_cut, -params.eta2)
                                      * fast::cos(params.k * r_cut - params.phi));
                pair_eng -= r_cut_eta1_arg + r_cut_eta2_arg;
                }

            return true;
            }
        else
            {
            return false;
            }
        }

    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("opp");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;        /// Stored rsq from the constructor
    Scalar radcontact; ///~ Stored contact-distance from the constructor [RHEOINF]
    unsigned int pair_typeids;///~ Stored pair typeIDs from the constructor [RHEOINF]
    unsigned int typei;///~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;///~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq;     /// Stored rcutsq from the constructor
    param_type params; /// Stored pair parameters for a given type pair
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_OPP_H__
