//~ ########## Created by the PRO-CF research group ##########
////~ HOOMD-blue:
//// Copyright (c) 2009-2022 The Regents of the University of Michigan.
//// Part of HOOMD-blue, released under the BSD 3-Clause License.
////~
////~ This file:
/////~ Written by Rob Campbell (2024) based on EvaluatorPairMorse.h


// ########## Modified by PRO-CF //~ [PROCF2023] [PROCF2024] ##########

#ifndef __PAIR_EVALUATOR_MORSE_H__
#define __PAIR_EVALUATOR_MORSE_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairMorseRepulse.h
    \brief Defines the pair evaluator class for Morse potential
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

    EvaluatorPairMorseRepulse evaluates the sum of two functions, the Morse potential:
    \f[ V_{\mathrm{Morse}}(r) = D_0 \left[ \exp \left(-2\alpha \left(r - r_0\right) \right)
                                           -2\exp \left(-\alpha \left(r-r_0\right) \right)  \right]
    \f]

    and one of two different types of repulsion: electrostatic repulsion (a la DLVO), or Yukawa repulsion:
    \f[ V_{\mathrm{electrostatic}}(r) = Z \frac{ (a_i a_j)}{r_0} \exp \left( -kappa_e \left(r - r_0\right) \right)
    \f]

    \f[ V_{\mathrm{yukawa}}(r) = \varepsilon \frac{ \exp \left( -\kappa_y r \right) }{r} \f]    

    NOTE: The added repulsion will decrease the depth of the attractive potential well. If you want to achieve the 
          same attraction as in a pure Morse case with the addition of repulsion, then you must adjust the parameters
          to be sure the potential energy curve has the correct form
*/

class EvaluatorPairMorseRepulse //~ change name to MorseRepulse
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        bool Erep; //~ add bool for Electrostatics on/off [PROCF2024]
        Scalar Z; //~ add Z param [PROCF2024]
        Scalar kappa_e; //~ add kappa_e param [PROCF2024]
        bool Yrep; //~ add bool for Yukawa repulsion on/off [PROCF2024]
        Scalar epsilon; //~ add epsilon param [PROCF2024]
        Scalar kappa_y; //~ add kappa_y param [PROCF2024]
        Scalar a_i; //~ add a_i param [PROCF2024]
        Scalar a_j; //~ replace r0 with a_j param [PROCF2024]
        Scalar f_contact; //~ add f_contact param [PROCF2023]


        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // CUDA memory hints
        void set_memory_hints() const { }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), Erep(false), Z(0), kappa_e(0), Yrep(false), epsilon(0), kappa_y(0), a_i(0), a_j(0), f_contact(0) { } //~ add params [PROCF2024]

        param_type(pybind11::dict v, bool managed = false, bool Erep = false, bool Yrep = false)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            this->Erep = Erep; //~ add boolean check for Electrostatic repulsion [PROCF2024]
            Z = v["Z"].cast<Scalar>(); //~ add Z param [PROCF2024]
            kappa_e = v["kappa_e"].cast<Scalar>(); //~ add kappa_e param [PROCF2024]
            this->Yrep = Yrep; //~ add boolean check for Yukawa repulsion [PROCF2024]
            epsilon = v["epsilon"].cast<Scalar>(); //~ add epsilon param [PROCF2024]
            kappa_y = v["kappa_y"].cast<Scalar>(); //~ add kappa_y param [PROCF2024]
            a_i = v["a_i"].cast<Scalar>(); //~ add a_i param [PROCF2024]
            a_j = v["a_j"].cast<Scalar>(); //~ add replace r0 with a_j param [PROCF2024]
            f_contact = v["f_contact"].cast<Scalar>(); //~ add f_contact param [PROCF2023]
            }

        param_type(Scalar d, Scalar a, bool Erep, Scalar z, Scalar ke, bool Yrep, Scalar e, Scalar ky, Scalar ai, Scalar aj, Scalar f, bool managed = false) //~ add params [PROCF2024]

            {
            D0 = d;
            alpha = a;
            Erep = Erep; //~ add boolean check for Electrostatic repulsion [PROCF2024]
            Z = z; //~ add Z param [PROCF2024]
            kappa_e = ke; //~ add kappa_e param [PROCF2024]
            Yrep = Yrep; //~ add boolean check for Yukawa repulsion [PROCF2024]
            epsilon = e; //~ add epsilon param [PROCF2024]
            kappa_y = ky; //~ add kappa_y param [PROCF2024]
            a_i = ai; //~ add a_i param [PROCF2024]
            a_j = aj; //~ replace r0 with a_j param [PROCF2024]
            f_contact = f; //~ add f_contact param [PROCF2023]
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["Erep"] = Erep; //~ add boolean check for Electrostatic repulsion [PROCF2024]
            v["Z"] = Z; //~ add Z param [PROCF2024]
            v["kappa_e"] = kappa_e; //~ add kappa_e param [PROCF2024]
            v["Yrep"] = Yrep; //~ add boolean check for Yukawa repulsion [PROCF2024]
            v["epsilon"] = epsilon; //~ add epsilon param [PROCF2024]
            v["kappa_y"] = kappa_y; //~ add kappa_y param [PROCF2024]
            v["a_i"] = a_i; //~ add a_i param [PROCF2024]
            v["a_j"] = a_j; //~ replace r0 with a_j param [PROCF2024]
            v["f_contact"] = f_contact; //~ add f_contact param [PROCF2023]
 
            return v;
            }
#endif
        } __attribute__((aligned(16)));

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairMorseRepulse(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), D0(_params.D0), alpha(_params.alpha), Erep(_params.Erep), Z(_params.Z), kappa_e(_params.kappa_e), Yrep(_params.Yrep), epsilon(_params.epsilon), kappa_y(_params.kappa_y), a_i(_params.a_i), a_j(_params.a_j), f_contact(_params.f_contact) //~ add params [PROCF2024]
        {
        }

    //! Morse doesn't use charge
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
        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar radsum = a_i + a_j; //~ r0 = radsum
            Scalar r = fast::sqrt(rsq);
            Scalar Exp_factor = fast::exp(-alpha * (r - radsum));

            //~ add contact force [PROCF2023]
            //~ if particles overlap (r < r0) apply contact force
            if(r < radsum)force_divr = f_contact * (Scalar(1.0) - (r-radsum)) * pow((Scalar(0.50)*radsum),3) / r;

            else{
                //~ calculate force as normal
                force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;

                //~ add repulsion
                if (Erep) {
                    // add electrostatic repulsion to the force divided by r in force_divr
                    Scalar rcutinv = fast::rsqrt(rcutsq);
                    Scalar rcut = Scalar(1.0) / rcutinv;
                    if (r < rcut && kappa_e != 0)
                        {
                        Scalar radprod = a_i * a_j; 
                        Scalar rmds = r - radsum;
                        Scalar radsuminv = Scalar(1.0) / radsum;

                        Scalar exp_val_e = fast::exp(-kappa_e * rmds);

                        Scalar forcerep_divr = kappa_e * radprod * radsuminv * Z * exp_val_e / r;
                        force_divr += forcerep_divr;
                        pair_eng += r * forcerep_divr / kappa_e;
                        }
                    }

                if (Yrep) {            
                  Scalar rinv = fast::rsqrt(rsq);
                  Scalar r2inv = Scalar(1.0) / rsq;
                  Scalar exp_val_y = fast::exp(-kappa_y * r);
            
                  force_divr += epsilon * exp_val_y * r2inv * (rinv + kappa_y);
                  pair_eng += epsilon * exp_val_y * rinv;
                  }

                //~ but still include contact force within 0.001 dist of colloid-colloid contact 
                Scalar Del_max = Scalar(0.001); //~ 0.001 or 0.01
                if(r<(radsum+Del_max))force_divr += f_contact * pow((Scalar(1.0) - (r-radsum)/Del_max), 3) * pow((Scalar(0.50)*radsum),3) / r;

                } 
            //~ 

            pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
            //~ force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r; //~ move this into overlap check [PROCF2023]

            //~ turn off energy shift [PROCF2024]
            //~ form must change because of new repulsion; but also shouldn't be needed because of contact force [PROCF2024]
            //if (energy_shift)
            //    {
            //    Scalar rcut = fast::sqrt(rcutsq);
            //    Scalar Exp_factor_cut = fast::exp(-alpha * (rcut - r0));
            //    pair_eng -= D0 * Exp_factor_cut * (Exp_factor_cut - Scalar(2.0));
            //    }
            //~
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
        return std::string("morse");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar D0;     //!< Depth of the Morse potential at its minimum
    Scalar alpha;  //!< Controls width of the potential well
    bool Erep;        //!<~ Electrostatic on/off bool [PROCF2024]
    Scalar Z;         //!<~ (scaled) surface charge [PROCF2024]
    Scalar kappa_e;   //!<~ Debye (screeneing) length [PROCF2024]
    bool Yrep;        //!<~ Yukawa on/off bool [PROCF2024]
    Scalar epsilon;   //!<~ Yukawa energy factor [PROCF2024]
    Scalar kappa_y;   //!<~ Yukawa scaling factor [PROCF2024]
    Scalar a_i;       //!<~ particle i radius [PROCF2024]
    Scalar a_j;       //!<~ particle j radius (replaces r0) [PROCF2024]
    //~Scalar r0;     //!< Offset, i.e., position of the potential minimum //~ comment out [PROCF2024]
    Scalar f_contact; //!< Contact force magnitude, for resolving overlap [PROCF2023]
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_H__
