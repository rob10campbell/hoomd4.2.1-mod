//~ ########## Created by the Rheoinformatic research group ##########
////~ HOOMD-blue:
//// Copyright (c) 2009-2022 The Regents of the University of Michigan.
//// Part of HOOMD-blue, released under the BSD 3-Clause License.
////~
////~ This file:
/////~ Written by Rob Campbell (2024) based on EvaluatorPairMorse.h


// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_MORSE_BRUSH2_H__
#define __PAIR_EVALUATOR_MORSE_BRUSH2_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairMorseBrush2.h
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

    EvaluatorPairMorseBrush2 evaluates the sum of two functions, the Morse potential:
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

class EvaluatorPairMorseBrush2 //~ change name to MorseBrush2
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        Scalar hb;
        Scalar soft_shift;
        Scalar Z; //~ add Z param [RHEOINF]
        Scalar kappa_e; //~ add kappa_e param [RHEOINF]
        Scalar f_contact; //~ add f_contact param [RHEOINF]
        bool scaled_D0; //~ add scaled_D0 param [RHEOINF]


        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // CUDA memory hints
        void set_memory_hints() const { }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), hb(0), soft_shift(0), Z(0), kappa_e(0), f_contact(0), scaled_D0(false) { } //~ add params [RHEOINF]

        param_type(pybind11::dict v, bool managed = false, bool Erep = false, bool Yrep = false)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            hb = v["hb"].cast<Scalar>();
            soft_shift = v["soft_shift"].cast<Scalar>();
            Z = v["Z"].cast<Scalar>(); //~ add Z param [RHEOINF]
            kappa_e = v["kappa_e"].cast<Scalar>(); //~ add kappa_e param [RHEOINF]
            f_contact = v["f_contact"].cast<Scalar>(); //~ add f_contact param [RHEOINF]
            this->scaled_D0 = scaled_D0; //~ add scaled_D0 param [RHEOINF
            }

        param_type(Scalar d, Scalar a, Scalar h, Scalar ss, Scalar z, Scalar ke, Scalar f, bool sD, bool managed = false) //~ add params [RHEOINF]

            {
            D0 = d;
            alpha = a;
            hb = h;
            soft_shift = ss;
            Z = z; //~ add Z param [RHEOINF]
            kappa_e = ke; //~ add kappa_e param [RHEOINF]
            f_contact = f; //~ add f_contact param [RHEOINF]
            scaled_D0 = sD; //~ add scaled_D0 param [RHEOINF]
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["hb"] = hb;
            v["soft_shift"] = soft_shift;
            v["Z"] = Z; //~ add Z param [RHEOINF]
            v["kappa_e"] = kappa_e; //~ add kappa_e param [RHEOINF]
            v["f_contact"] = f_contact; //~ add f_contact param [RHEOINF]
            v["scaled_D0"] = scaled_D0; //~ add scaled_D0 param [RHEOINF] 
            return v;
            }
#endif
        } __attribute__((aligned(16)));

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairMorseBrush2(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~ add radcontact and pair typeids [RHEOINF]
        : rsq(_rsq), rcutsq(_rcutsq), D0(_params.D0), alpha(_params.alpha), hb(_params.hb), soft_shift(_params.soft_shift), Z(_params.Z), kappa_e(_params.kappa_e), f_contact(_params.f_contact), scaled_D0(_params.scaled_D0) //~ add params [RHEOINF]
        {
        typei = _pair_typeids[0]; //~ add typei [RHEOINF]
        typej = _pair_typeids[1]; //~ add typej [RHEOINF]
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
        //~ Add radsum from passed diameters [RHEOINF] 
        Scalar radsum_brush = 0.5 * (diameter_i + diameter_j) + (2*hb); 
        Scalar radsum_colloid = 0.5 * (diameter_i + diameter_j); 
        //~ Scale attraction strength by particle size is scaled_D0 is true
        if (scaled_D0)
          {
          D0 = D0 * (0.5*radsum_brush);
          }   
        //~ 

        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq);

            //~ add repulsion
            if (kappa_e)
                {
                Scalar radprod = (Scalar(0.5)*diameter_i) * (Scalar(0.5)*diameter_j); 
                Scalar rmds = r - radsum_colloid;
                Scalar radsuminv = Scalar(1.0) / radsum_colloid;

                Scalar exp_val_e = fast::exp(-kappa_e * rmds);

                Scalar forcerep_divr = kappa_e * radprod * radsuminv * Z * exp_val_e / r;
                force_divr += forcerep_divr;
                pair_eng += r * forcerep_divr / kappa_e;
                }


            // add morse with brush
            if (r < radsum_brush)
                {
                Scalar alpha_shift = soft_shift * alpha;
                Scalar Exp_factor_shift = fast::exp(-alpha_shift * (r - radsum_brush));
                force_divr += Scalar(2.0) * D0 * alpha_shift * Exp_factor_shift * (Exp_factor_shift - Scalar(1.0)) / r;
                pair_eng += D0 * Exp_factor_shift * (Exp_factor_shift - Scalar(2.0));
                
                if (energy_shift)
                    {
                    Scalar rcut = fast::sqrt(rcutsq);
                    Scalar Exp_factor_shift_cut = fast::exp(-alpha_shift * (rcut - radsum_brush));
                    pair_eng -= D0 * Exp_factor_shift_cut * (Exp_factor_shift_cut - Scalar(2.0));
                    }
                }

            else
                {
                Scalar Exp_factor = fast::exp(-alpha * (r - radsum_brush));
                force_divr += Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;
                pair_eng += D0 * Exp_factor * (Exp_factor - Scalar(2.0));

                if (energy_shift)
                    {
                    Scalar rcut = fast::sqrt(rcutsq);
                    Scalar Exp_factor_cut = fast::exp(-alpha * (rcut - radsum_brush));
                    pair_eng -= D0 * Exp_factor_cut * (Exp_factor_cut - Scalar(2.0));
                    }

                }

            //std::cout << "--MorseBrush2--" << std::endl;
            //std::cout << "Z: " << Z << std::endl;
            //std::cout << "kappa_e: " << kappa_e << std::endl;
            //std::cout << "diameter_i: " << diameter_i << std::endl;
            //std::cout << "D0: " << D0 << std::endl;
            //std::cout << "hb: " << hb << std::endl;
            //std::cout << "pair_eng: " << pair_eng << std::endl;
            //std::cout << "force_divr: " << force_divr << std::endl;
            //std::cout << "-------------" << std::endl;

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
        return std::string("morsebrush2");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar radcontact;//!< Stored contact-distance from the constructor [RHEOINF]
    unsigned int pair_typeids;//!< Stored pair typeIDs from the constructor [RHEOINF]
    unsigned int typei;//!<~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;//!<~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar diameter_i;//!<~ add diameter_i [RHEOINF]
    Scalar diameter_j;//!<~ add diameter_j [RHEOINF]
    Scalar D0;     //!< Depth of the Morse potential at its minimum
    Scalar alpha;  //!< Controls width of the potential well
    Scalar hb; 
    Scalar soft_shift;
    Scalar Z;         //!<~ (scaled) surface charge [RHEOINF]
    Scalar kappa_e;   //!<~ Debye (screeneing) length [RHEOINF]
    Scalar f_contact; //!< Contact force magnitude, for resolving overlap [RHEOINF]
    bool scaled_D0;   //!<~ on/off bool for scaling D0 by particle size [RHEOINF]
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_BRUSH2_H__
