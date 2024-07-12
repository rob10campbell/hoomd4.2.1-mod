// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_REACTION_FIELD_H__
#define __PAIR_EVALUATOR_REACTION_FIELD_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h" // add vectors for optional position [RHEOINF]

/*! \file EvaluatorPairReactionField.h
    \brief Defines the pair evaluator class for ReactionField potentials
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
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
//! Class for evaluating the Onsager reaction field pair potential
/*! <b>General Overview</b>

    See EvaluatorPairLJ

    <b>ReactionField specifics</b>

    EvaluatorPairReactionField evaluates the function:
    \f[ V_{\mathrm{RF}}(r) = \varepsilon \left[ \frac{1}{r} +
        \frac{(\epsilon_{RF}-1) r^2}{(2 \epsilon_{RF} + 1) r_c^3} \right]\f]

    If \epsilon_{RF} is zero, it will be treated as infinity.
*/
class EvaluatorPairReactionField
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        // potential parameters
        Scalar eps, eps_rf;
        bool use_charge;

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // set CUDA memory hints
        void set_memory_hint() const { }
#endif

#ifndef __HIPCC__
        param_type() : eps(0), eps_rf(0), use_charge(false) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            eps = v["epsilon"].cast<Scalar>();
            eps_rf = v["eps_rf"].cast<Scalar>();
            use_charge = v["use_charge"].cast<bool>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["epsilon"] = eps;
            v["eps_rf"] = eps_rf;
            v["use_charge"] = use_charge;

            return v;
            }
#endif
        } __attribute((aligned(16)));

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _radcontact the sum of the interacting particle radii [RHEOINF]
        \param _pair_typeids the typeIDs of the interacting particles [RHEOINF]
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairReactionField(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~ add radcontact, pair_typeIDs [RHEOINF]
        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), epsilon(_params.eps), epsrf(_params.eps_rf), //~ add radcontact [RHEOINF] 
          use_charge(_params.use_charge), qiqj(1.0)
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

    //! ReactionField uses charge
    DEVICE static bool needsCharge()
        {
        return true;
        }

    //! Accept the optional charge values
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj)
        {
        if (use_charge)
            qiqj = qi * qj;
        }

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
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that V(r) is continuous at the
       cutoff \note There is no need to check if rsq < rcutsq in this method. Cutoff tests are
       performed in PotentialPair.

        \return True if they are evaluated or false if they are not because we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq && epsilon != 0 && qiqj != 0)
            {
            Scalar rcut3inv = fast::rsqrt(rcutsq) / rcutsq;
            Scalar rinv = fast::rsqrt(rsq);
            Scalar r = Scalar(1.0) / rinv;
            Scalar r2inv = Scalar(1.0) / rsq;

            Scalar eps_fac = (epsrf - Scalar(1.0)) / (Scalar(2.0) * epsrf + Scalar(1.0)) * rcut3inv;
            if (epsrf == Scalar(0.0))
                {
                eps_fac = Scalar(1.0 / 2.0) * rcut3inv;
                }

            force_divr = qiqj * epsilon * (r2inv * rinv - Scalar(2.0) * eps_fac);
            pair_eng = qiqj * epsilon * (rinv + eps_fac * r * r);

            if (energy_shift)
                {
                Scalar rcutinv = fast::rsqrt(rcutsq);
                Scalar rcut = Scalar(1.0) / rcutinv;
                pair_eng -= qiqj * epsilon * (rcutinv + eps_fac * rcut * rcut);
                }
            return true;
            }
        else
            return false;
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
        return std::string("reaction_field");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;      //!< Stored rsq from the constructor
    Scalar radcontact;//!< Stored contact-distance from the constructor [RHEOINF]
    unsigned int pair_typeids;//!< Stored pair typeIDs from the constructor [RHEOINF]
    unsigned int typei;//!<~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;//!<~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq;   //!< Stored rcutsq from the constructor
    Scalar epsilon;  //!< epsilon parameter extracted from the params passed to the constructor
    Scalar epsrf;    //!< epsilon_rf parameter extracted from the params passed to the constructor
    bool use_charge; //!< True if we are using the particle charges
    Scalar qiqj;     //!< Product of charges
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_REACTION_FIELD_H__
