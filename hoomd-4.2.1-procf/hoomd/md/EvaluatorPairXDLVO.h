//~ ########## Created by the Rheoinformatic research group ##########
////~ HOOMD-blue:
//// Copyright (c) 2009-2022 The Regents of the University of Michigan.
//// Part of HOOMD-blue, released under the BSD 3-Clause License.
////~
////~ This file:
/////~ Written by Rob Campbell (2024) based on EvaluatorPairMorse.h


// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_XDLVO_H__
#define __PAIR_EVALUATOR_XDLVO_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairXDLVO.h
    \brief Defines the pair evaluator class for xDLVO potential
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
//! Class for evaluating the xDLVO pair potential
/*! <b>General Overview</b>

    See EvaluatorPairLJ.

    <b>xDLVO specifics</b>

    EvaluatorPairXDLVO evaluates the sum of two functions, the Morse potential:
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

class EvaluatorPairXDLVO //~ change name to XDLVO 
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        Scalar Z; //~ add Z param [RHEOINF]
        Scalar kappa_e; //~ add kappa_e param [RHEOINF]
        Scalar B; //~ add B param [RHEOINF]
        Scalar r_e; //~ add r_e param [RHEOINF]
        Scalar rho_e; //~ add rho_e param [RHEOINF]
        Scalar sigma; //~ add sigma param [RHEOINF]
        Scalar kT; //~ add kT param [RHEOINF]
        Scalar f_contact; //~ add f_contact param [RHEOINF]
        bool scaled_D0; //~ add scaled_D0 param [RHEOINF]


        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // CUDA memory hints
        void set_memory_hints() const { }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), Z(0), kappa_e(0), B(0), r_e(0), rho_e(0), sigma(0), kT(0), f_contact(0), scaled_D0(false) { } //~ add params [RHEOINF]

        param_type(pybind11::dict v, bool managed = false, bool Erep = false, bool Yrep = false)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            Z = v["Z"].cast<Scalar>(); //~ add Z param [RHEOINF]
            kappa_e = v["kappa_e"].cast<Scalar>(); //~ add kappa_e param [RHEOINF]
            B = v["B"].cast<Scalar>(); //~ add B param [RHEOINF]
            r_e = v["r_e"].cast<Scalar>(); //~ add r_e param [RHEOINF]
            rho_e = v["rho_e"].cast<Scalar>(); //~ add rho_e param [RHEOINF]
            sigma = v["sigma"].cast<Scalar>(); //~ add sigma param [RHEOINF]
            kT = v["kT"].cast<Scalar>(); //~ add kT param [RHEOINF]
            f_contact = v["f_contact"].cast<Scalar>(); //~ add f_contact param [RHEOINF]
            this->scaled_D0 = scaled_D0; //~ add scaled_D0 param [RHEOINF
            }

        param_type(Scalar d, Scalar a, Scalar z, Scalar ke, Scalar b, Scalar re, Scalar rh, Scalar sg, Scalar kt, Scalar f, bool sD, bool managed = false) //~ add params [RHEOINF]

            {
            D0 = d;
            alpha = a;
            Z = z; //~ add Z param [RHEOINF]
            kappa_e = ke; //~ add kappa_e param [RHEOINF]
            B = b; //~ add B param [RHEOINF]
            r_e = re; //~ add r_e param [RHEOINF]
            rho_e = rh; //~ add rho_e param [RHEOINF]
            sigma = sg; //~ add sigma param [RHEOINF]
            kT = kt; //~ add kT param [RHEOINF]
            f_contact = f; //~ add f_contact param [RHEOINF]
            scaled_D0 = sD; //~ add scaled_D0 param [RHEOINF]
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["Z"] = Z; //~ add Z param [RHEOINF]
            v["kappa_e"] = kappa_e; //~ add kappa_e param [RHEOINF]
            v["B"] = B; //~ add B param [RHEOINF]
            v["r_e"] = r_e; //~ add r_e param [RHEOINF]
            v["rho_e"] = rho_e; //~ add rho_e param [RHEOINF]
            v["sigma"] = sigma; //~ add sigma param [RHEOINF]
            v["kT"] = kT; //~ add kT param [RHEOINF]
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
    DEVICE EvaluatorPairXDLVO(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~ add radcontact and pair typeids [RHEOINF]
        : rsq(_rsq), rcutsq(_rcutsq), D0(_params.D0), alpha(_params.alpha), Z(_params.Z), kappa_e(_params.kappa_e), B(_params.B), r_e(_params.r_e), rho_e(_params.rho_e), sigma(_params.sigma), kT(_params.kT), f_contact(_params.f_contact), scaled_D0(_params.scaled_D0) //~ add params [RHEOINF]
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
        Scalar radsum = 0.5 * (diameter_i + diameter_j); 
        //~ Scale attraction strength by particle size is scaled_D0 is true
        if (scaled_D0)
          {
          D0 = D0 * (0.5*radsum);
          }   
        //~ 

        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq);
            Scalar Exp_factor = fast::exp(-alpha * (r - radsum));

            //~ check if contact force is provided [RHEOINF]
            if (f_contact != 0.0)
                {
                //~ if particles overlap (r < radsum) apply contact force
                if(r < radsum)force_divr = f_contact * (Scalar(1.0) - (r-radsum)) * pow((Scalar(0.50)*radsum),3) / r;

                //~ if particles do not overlap
                else{
                    //~ calculate Morse force as normal
                    force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;

                    //~ but still include contact force within 0.001 dist of colloid-colloid contact 
                    Scalar Del_max = Scalar(0.001); //~ 0.001 or 0.01
                    if(r<(radsum+Del_max))force_divr += f_contact * pow((Scalar(1.0) - (r-radsum)/Del_max), 3) * pow((Scalar(0.50)*radsum),3) / r;
                    } 

                }

            //if f_contact=0.0
            else{
                //~ calculate force as normal (use Morse repulsion)
                force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;
                }

            pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));


            //~ apply repulsion WITHIN A PORTION of the Morse range
            //Scalar rcut = fast::sqrt(rcutsq);
            //if (r < Scalar(0.75)*rcut) 
            //   {

            // add electrostatic repulsion to the force divided by r in force_divr
            //Scalar rcutinv = fast::rsqrt(rcutsq);
            //Scalar rcut = Scalar(1.0) / rcutinv;
            //if (r < rcut && kappa_e != 0)
            if (kappa_e != 0)
                {
                Scalar radprod = (Scalar(0.5)*diameter_i) * (Scalar(0.5)*diameter_j); 
                Scalar rmds = r - radsum;
                Scalar radsuminv = Scalar(1.0) / radsum;

                Scalar exp_val_e = fast::exp(-kappa_e * rmds);

                Scalar forcerep_divr = kappa_e * radprod * radsuminv * Z * exp_val_e / r;
                force_divr += forcerep_divr;
                pair_eng += r * forcerep_divr / kappa_e;
                }



            // add ionic osmotic potential for salt effects (Donnan potential)
            Scalar D_ij = radsum + 2*r_e + 2*sigma; //~ center-center 2-particle range of ionic effects
            Scalar d_hydroshell = radsum + 2*sigma; //~ particle + water/hydration layer
            if (r >= d_hydroshell && r <= 2*D_ij)
                {
                Scalar forceosm_divr = -Scalar(0.25) * M_PI * kT * rho_e * (4*pow(D_ij,2) - pow(r,2)) / r;
                force_divr += forceosm_divr;
                pair_eng += -Scalar(4/3) * M_PI * kT * pow(D_ij,3) * rho_e * ( 1 - ((3*r)/(4*D_ij)) + (pow(r,3)/(16*pow(D_ij,3))) );
                }




            // add ion-particle interaction potential for counter-ion bridging
            if ((r - radsum) <= 0.006)
              {
              // assuming a 1-1 salt like NaCl
              Scalar z_ion = 1;
              Scalar ion_radprod = (Scalar(0.5)*diameter_i) * r_e;
              Scalar ion_radsum = Scalar(0.5)*diameter_i + r_e;
              Scalar ion_radsuminv = Scalar(1.0) / ion_radsum;

              // PARTICLE I
              Scalar hi_min = Scalar(0.5)*diameter_i + sigma;
              Scalar hi_max = hi_min + (r - radsum); //r_e;

              // Numerical integration using the trapezoidal rule
              Scalar integral_i = 0.0;
              Scalar stepsize_i = r_e; //0.0001; 
              Scalar steps_i = std::round((hi_max - hi_min) / stepsize_i);
              //Scalar stepsize_i = (hi_max - hi_min) / steps_i; // Step size

              // intermediate points : step_size*Sum(f(x_i))
              // force is only calculated to find energy at each point
              for (int i = 1; i < steps_i; ++i)
                {
                Scalar r_ionshell = hi_min + (i * stepsize_i);
                Scalar rmds_i = r_ionshell - (Scalar(0.5)*diameter_i + r_e);
                Scalar exp_val_e_i = fast::exp(-kappa_e * rmds_i);
                Scalar force_i = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_i;
                Scalar electro_eng_i = force_i / kappa_e;
                integral_i += (rho_e * std::exp(-z_ion * electro_eng_i / kT)) / r_ionshell;
                }
              integral_i *= stepsize_i;

              // min boundary
              // force is only calculated to find energy
              Scalar rmds_min = hi_min - ion_radsum;
              Scalar exp_val_e_min = fast::exp(-kappa_e * rmds_min);
              Scalar force_min = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_min;
              Scalar electro_eng_min = force_min / kappa_e;
              Scalar min_i = (rho_e * std::exp(-z_ion * electro_eng_min / kT) / hi_min);

              // max boundary
              // force is only calculated to find energy
              Scalar rmds_max = hi_max - ion_radsum;
              Scalar exp_val_e_max = fast::exp(-kappa_e * rmds_max);
              Scalar force_max = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_max;
              Scalar electro_eng_max = force_max / kappa_e;
              Scalar max_i = (rho_e * std::exp(-z_ion * electro_eng_max / kT) / hi_max);

              // boundary points : (step_size/2) * (f(x_0)+f(x_N))
              integral_i += Scalar(0.5) * stepsize_i * (min_i+max_i);

              Scalar totali_eng = Scalar(4.0) * M_PI * B * integral_i;

              // use central difference formula to numerically approximate dU/dr = force
              Scalar deltaR_force = 0.0001;

              Scalar hR_min = hi_min + deltaR_force;
              Scalar rmds_Rmin = hR_min - ion_radsum;
              Scalar exp_val_e_Rmin = fast::exp(-kappa_e * rmds_Rmin);
              Scalar force_Rmin = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_Rmin;
              Scalar Ue_plus = force_Rmin / kappa_e;

              Scalar hR_max = hi_max - deltaR_force;
              Scalar rmds_Rmax = hR_max - ion_radsum;
              Scalar exp_val_e_Rmax = fast::exp(-kappa_e * rmds_Rmax);
              Scalar force_Rmax = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_Rmax;
              Scalar Ue_minus = force_Rmax / kappa_e;

              Scalar force_ioni_divr = -(Ue_plus - Ue_minus) / (2*deltaR_force);

              // PARTICLE J
              Scalar hj_min = Scalar(0.5)*diameter_j + sigma;
              Scalar hj_max = hi_min + (r - radsum); //r_e;

              // Numerical integration using the trapezoidal rule
              Scalar integral_j = 0.0;
              Scalar stepsize_j =r_e; //0.0001; 
              Scalar steps_j = (hj_max - hj_min) / stepsize_j;
              //Scalar stepsize_j = (hj_max - hj_min) / steps_j; // Step size

              // intermediate points : step_size*Sum(f(x_i))
              // force is only calculated to find energy at each point
              for (int i = 1; i < steps_j; ++i)
                {
                Scalar r_ionshell = hj_min + (i * stepsize_j);
                Scalar rmds_j = r_ionshell - (Scalar(0.5)*diameter_j + r_e);
                Scalar exp_val_e_j = fast::exp(-kappa_e * rmds_j);
                Scalar force_j = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_j;
                Scalar electro_eng_j = force_j / kappa_e;
                integral_j += (rho_e * std::exp(-z_ion * electro_eng_j / kT)) / r_ionshell;
                }
              integral_j *= stepsize_j;

              // min boundary
              // force is only calculated to find energy
              Scalar rmdsj_min = hj_min - ion_radsum;
              Scalar exp_val_e_jmin = fast::exp(-kappa_e * rmdsj_min);
              Scalar force_jmin = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_jmin;
              Scalar electro_eng_jmin = force_jmin / kappa_e;
              Scalar min_j = (rho_e * std::exp(-z_ion * electro_eng_jmin / kT) / hj_min);

              // max boundary
              // force is only calculated to find energy
              Scalar rmdsj_max = hj_max - ion_radsum;
              Scalar exp_val_e_jmax = fast::exp(-kappa_e * rmdsj_max);
              Scalar force_jmax = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_jmax;
              Scalar electro_eng_jmax = force_jmax / kappa_e;
              Scalar max_j = (rho_e * std::exp(-z_ion * electro_eng_jmax / kT) / hj_max);

              // boundary points : (step_size/2) * (f(x_0)+f(x_N))
              integral_j += Scalar(0.5) * stepsize_j * (min_j+max_j);

              Scalar totalj_eng = Scalar(4.0) * M_PI * B * integral_j;

              // use central difference formula to numerically approximate dU/dr = force
              //Scalar deltaR_force = 0.0001;

              Scalar hRj_min = hj_min + deltaR_force;
              Scalar rmds_Rjmin = hRj_min - ion_radsum;
              Scalar exp_val_e_Rjmin = fast::exp(-kappa_e * rmds_Rjmin);
              Scalar force_Rjmin = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_Rjmin;
              Scalar Uej_plus = force_Rjmin / kappa_e;

              Scalar hRj_max = hj_max - deltaR_force;
              Scalar rmds_Rjmax = hRj_max - ion_radsum;
              Scalar exp_val_e_Rjmax = fast::exp(-kappa_e * rmds_Rjmax);
              Scalar force_Rjmax = kappa_e * ion_radprod * ion_radsuminv * Z * exp_val_e_Rjmax;
              Scalar Uej_minus = force_Rjmax / kappa_e;

              Scalar force_ionj_divr = -(Uej_plus - Uej_minus) / (2*deltaR_force);

              force_divr += (force_ioni_divr - force_ionj_divr);
              pair_eng -= (totali_eng + totalj_eng);
              }
            
            //} // close second dist check

            //pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
            //~ force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r; //~ move this into overlap check [RHEOINF]

            //~ turn off energy shift [RHEOINF]
            //~ form must change because of new repulsion; but also shouldn't be needed because of contact force [RHEOINF]
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
        return std::string("xdlvo");
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
    Scalar Z;         //!<~ (scaled) surface charge [RHEOINF]
    Scalar kappa_e;   //!<~ Debye (screeneing) length [RHEOINF]
    Scalar B;   //!<~ particle-ion dispersion force [RHEOINF]
    Scalar r_e;     //!< hydrated ion size [RHEOINF]
    Scalar rho_e;   //!<~ charge/ion/salt density/concentration [RHEOINF]
    Scalar sigma;     //!< size of the hydration layer on the particle [RHEOINF]
    Scalar kT;     //!< system temperature [RHEOINF]
    Scalar f_contact; //!< Contact force magnitude, for resolving overlap [RHEOINF]
    bool scaled_D0;   //!<~ on/off bool for scaling D0 by particle size [RHEOINF]
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_XDLVO_H__
