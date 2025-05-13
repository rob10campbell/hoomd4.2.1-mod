//~ ########## Created by the Rheoinformatic research group ##########
////~ HOOMD-blue:
//// Copyright (c) 2009-2022 The Regents of the University of Michigan.
//// Part of HOOMD-blue, released under the BSD 3-Clause License.
////~
////~ This file:
/////~ Written by Rob Campbell (2025)
//

#ifndef __PAIR_EVALUATOR_MORSE_BRUSH_H__
#define __PAIR_EVALUATOR_MORSE_BRUSH_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairMorseBrush.h
    \brief Defines the pair evaluator class for Morse potential which has been modified to mimic the effect of a polymer brush coating as observed in experiments by Calvin Zhuang and Ali Mohraz at UC Irvine.
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
//! Class for evaluating the Morse pair potential
/*! <b>General Overview</b>

    See EvaluatorPairLJ.

    <b>Morse specifics</b>

    EvaluatorPairMorseBrush evaluates the function:
    \f[ V_{\mathrm{Morse}}(r) = D_0 \left[ \exp \left(-2\alpha \left(r - r_0\right) \right)
                                           -2\exp \left(-\alpha \left(r-r_0\right) \right)  \right]
   \f]

   If (r-r0)<0 then the evaluator applies a soft_shift to the repulsive part of the potential, alpha->soft_shift*alpha
   This allows the user to optionally soften the repulsive force, mimicking a change in density or rigidity of a polymer brush coating.

*/
class EvaluatorPairMorseBrush
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        Scalar hb;         //~ length of the polymer brush coating the colloid [RHEOINF]
        Scalar soft_shift; //~ the amount of shift in the soft repulsion (mimics density/stiffness of the brush) [RHEOINF]
        bool scaled_D0;    //~ add scaled_D0 param [RHEOINF]


        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // CUDA memory hints
        void set_memory_hints() const { }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), hb(0), soft_shift(0), scaled_D0(false) { } //~ add params [RHEOINF]

        param_type(pybind11::dict v, bool managed = false)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            hb = v["hb"].cast<Scalar>(); //~ [RHEOINF]
            soft_shift = v["soft_shift"].cast<Scalar>(); //~ [RHEOINF]
            this->scaled_D0 = scaled_D0; //~ [RHEOINF]
            }

        param_type(Scalar d, Scalar a, Scalar h, Scalar ss, bool sD, bool managed = false) //~ add params [RHEOINF]
            {
            D0 = d;
            alpha = a;
            hb = h; //~ [RHEOINF]
            soft_shift = ss; //~ [RHEOINF]
            scaled_D0 = sD; //~ [RHEOINF]
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["hb"] = hb; //~ [RHEOINF]
            v["soft_shift"] = soft_shift; //~ [RHEOINF]
            v["scaled_D0"] = scaled_D0; //~ [RHEOINF] 
            return v;
            }
#endif
        } __attribute__((aligned(16)));

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _radcontact the sum of the interacting particle radii [RHEOINF]
        \param _pair_typeids the typeIDs of the interacting particles [RHEOINF]
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairMorseBrush(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~add radcontact, pair_typeIDs [RHEOINF]
        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha),
         hb(_params.hb), soft_shift(_params.soft_shift), scaled_D0(_params.scaled_D0) //~ add radcontact, diameters, new params [RHEOINF]
        {
        }

    //~ add diameter [RHEOINF] 
    DEVICE static bool needsDiameter()
        {
        return true;
        }
    //! Accept the optional diameter values
    /*! \param di Diameter of particle i
        \param dj Diameter of particle j
    */
    DEVICE void setDiameter(Scalar di, Scalar dj)
        {
        diameter_i = di;
        diameter_j = dj;
        }
    //~

    //! MorseBrush doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional charge values.
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

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
        //~ Add radsum from passed diameters and include brush size for both particles [RHEOINF] 
        Scalar radsum_morse = 0.5 * (diameter_i + diameter_j) + (2*hb); 
        //Scalar rcut = fast::sqrt(rcutsq);

        //~ Scale attraction strength by particle size is scaled_D0 is true
        //~ NOTE: changing the value of a parameter passed from Python is bad form, but it does seem to work...
        if (scaled_D0)
          {
          D0 = D0 * (0.5*radsum_morse);
          }   

        // if particles are within range, compute the force divided by r and potential energy
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq); 

            //~ if particles overlap (r < radsum) use the SHIFTED Morse potential repulsion 
            if (r < radsum_morse)
                {
                Scalar alpha_shift = soft_shift * alpha;
                Scalar Exp_factor_shift = fast::exp(-alpha_shift * (r - radsum_morse));
                force_divr = Scalar(2.0) * D0 * alpha_shift * Exp_factor_shift * (Exp_factor_shift - Scalar(1.0)) / r;
                pair_eng = D0 * Exp_factor_shift * (Exp_factor_shift - Scalar(2.0));

                if (energy_shift)
                    {
                    Scalar rcut = fast::sqrt(rcutsq);
                    Scalar Exp_factor_shift_cut = fast::exp(-alpha_shift * (rcut - radsum_morse));
                    pair_eng -= D0 * Exp_factor_shift_cut * (Exp_factor_shift_cut - Scalar(2.0));
                    }
                }

            // otherwise calculate force as normal
            else
                {
                Scalar Exp_factor = fast::exp(-alpha * (r - radsum_morse));
                force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;
                pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));

                if (energy_shift)
                    {
                    Scalar rcut = fast::sqrt(rcutsq);
                    Scalar Exp_factor_cut = fast::exp(-alpha * (rcut - radsum_morse));
                    pair_eng -= D0 * Exp_factor_cut * (Exp_factor_cut - Scalar(2.0));
                    }

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
        return std::string("morsebrush");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;                 //!< Stored rsq from the constructor
    Scalar radcontact;          //!< Stored contact-distance from the constructor [RHEOINF]
    unsigned int pair_typeids;  //!< Stored pair typeIDs from the constructor [RHEOINF]
    unsigned int typei;         //!<~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;         //!<~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq;              //!< Stored rcutsq from the constructor
    Scalar diameter_i;          //!<~ add diameter_i [RHEOINF]
    Scalar diameter_j;          //!<~ add diameter_j [RHEOINF]
    Scalar D0;                  //!< Depth of the Morse potential at its minimum
    Scalar alpha;               //!< Controls width of the potential well
    Scalar hb;                  //!< Offset, i.e., position of the potential minimum
    Scalar soft_shift;          //!< Offset, i.e., position of the potential minimum
    bool scaled_D0;             //!<~ on/off bool for scaling D0 by particle size [RHEOINF]
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_BRUSH_H__
