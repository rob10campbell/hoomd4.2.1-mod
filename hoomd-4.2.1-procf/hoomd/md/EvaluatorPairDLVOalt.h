// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_DLVO_ALT_H__
#define __PAIR_EVALUATOR_DLVO_ALT_H__

#ifndef __HIPCC__
#include <string>
#include <cmath>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairDLVOalt.h
    \brief Defines the pair evaluator class for DLVO potential
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
//! DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
//! compiler
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
//! Class for evaluating the alternate electrostatic DLVO pair potential
/** See Israelachvili 2011, pp. 317.
 */
class EvaluatorPairDLVOalt
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar cs_mM;
        Scalar hb_nm;
        Scalar Zeta_mV;
        Scalar a1;
        Scalar a2;

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
        param_type() : cs_mM(0), hb_nm(0), Zeta_mV(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            cs_mM = v["cs_mM"].cast<Scalar>();
            hb_nm = v["hb_nm"].cast<Scalar>();
            Zeta_mV = v["Zeta_mV"].cast<Scalar>();
            a1 = v["a1"].cast<Scalar>();
            a2 = v["a2"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["cs_mM"] = cs_mM;
            v["hb_nm"] = hb_nm;
            v["Zeta_mV"] = Zeta_mV;
            v["a1"] = a1;
            v["a2"] = a2;
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
    DEVICE EvaluatorPairDLVOalt(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~ add radcontact, pair_typeIDs [RHEOINF]
        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), cs_mM(_params.cs_mM), hb_nm(_params.hb_nm), Zeta_mV(_params.Zeta_mV) //~ add radcontact [RHEOINF]
        {
        radsum = _params.a1 + _params.a2;
        radsub = _params.a1 - _params.a2;
        radprod = _params.a1 * _params.a2;
        radsumsq = _params.a1 * _params.a1 + _params.a2 * _params.a2;
        radsubsq = _params.a1 * _params.a1 - _params.a2 * _params.a2;
        delta = radsum - Scalar(1.0);
        typei = _pair_typeids[0]; //~ add typei [RHEOINF]
        typej = _pair_typeids[1]; //~ add typej [RHEOINF]
        }

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

    //! DLVO doesn't use charge
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
        // precompute some quantities
        Scalar rinv = fast::rsqrt(rsq);
        Scalar r = Scalar(1.0) / rinv;
        Scalar rcutinv = fast::rsqrt(rcutsq);
        Scalar rcut = Scalar(1.0) / rcutinv;

        // compute the force divided by r in force_divr
        if (r < rcut)
            {
            //Scalar sim_R_C = 1; // simulation units for particle size
            Scalar sim_kT = 0.1; // simulation units for kT
            //Scalar sim_energy = 4.1153799999999995e-20; // (k_B*T)/sim_kT

            // parameters in real units
            Scalar cs_Mol = cs_mM * 1e-3; // 10e-3
            Scalar hb = hb_nm * 1e-9; //12e-9
            Scalar Zeta = Zeta_mV * 1e-3; //50e-3

            // constants in real units
            Scalar epsilon_0 = 8.854e-12; // permittivity of free space [F/m]
            Scalar epsilon_r = 51.4;      // dielectric constant of 75w/w% glycerol water
            Scalar R_small = 0.5*1.1e-6;
            Scalar T = 298;
            Scalar k_B = 1.381e-23;
            Scalar e_charge = 1.602e-19;  // Coulombs (electron charge) [C]
            Scalar z1 = 1;
            Scalar Na = 6.022e23;
            Scalar rho_inf = 1000 * Na * cs_Mol;
            Scalar numerator = rho_inf * pow(e_charge,2) * z1 * 2;
            Scalar denomenator = epsilon_r * epsilon_0 * k_B * T;
            Scalar kappa = fast::rsqrt(numerator/denomenator);
            Scalar sigma = kappa * epsilon_0 * epsilon_r* Zeta;
            Scalar Z2R_over_kT = 2 * R_small * M_PI * pow(sigma,2)/(pow(kappa,3) * epsilon_0 * epsilon_r * (2*hb) * k_B*T);
            //Scalar Z2R_over_kT = 4 * (R_small/2) * M_PI * pow(sigma,2)/(pow(kappa,3) * epsilon_0 * epsilon_r * (2*hb) * k_B*T);
            //Scalar Z2R_over_kT = 4 * ((R_small*R_small)/(R_small+R_small)) * M_PI * pow(sigma,2)/(pow(kappa,3) * epsilon_0 * epsilon_r * (2*hb) * k_B*T);
            //Scalar Z_over_kT = 4 * M_PI * pow(sigma,2) / (pow(kappa,3) * epsilon_0 * epsilon_r * (2*hb) * k_B*T);
            //Scalar Z_over_kT_sim = Z_over_kT * sim_kT;
            Scalar Z2R_over_kT_sim = Z2R_over_kT * sim_kT;

            //Scalar radsum = (0.5*diameter_i)+(0.5*diameter_j);
            Scalar radsum = 2;
            //Scalar radprod = 1;
            //Scalar radsuminv = 0.5;
            //Scalar radprod = (0.5*diameter_i) * (0.5*diameter_j);
            Scalar h_ij = r - radsum;
            Scalar hb_sim = hb / R_small;

            if (h_ij < (2*hb_sim))
              {
              //force_divr = kappa * radprod * radsuminv * Z_sim * fast::exp(-kappa*h_ij) / r;
              force_divr = kappa * Z2R_over_kT_sim * fast::exp(-kappa*h_ij) / r;
              pair_eng = radsum * Z2R_over_kT_sim * (fast::exp(kappa*h_ij)-1) * fast::exp(-kappa*h_ij);
              }
            else
              {
              force_divr = kappa * Z2R_over_kT_sim * (fast::exp(kappa*(2*hb_sim))-1) * (-1*fast::exp(-kappa*h_ij)) / r;
              pair_eng = radsum * Z2R_over_kT_sim * (fast::exp(kappa*(2*hb_sim))-1) * fast::exp(-kappa*h_ij);
              }

            /*
            if (energy_shift)
                {
                Scalar rcutt = rcut;
                Scalar rmdscut = rcutt - radsum;
                Scalar rmdsqscut = rcutt * rcutt - radsum * radsum;
                Scalar rmdsqmcut = rcutt * rcutt - radsub * radsub;
                Scalar rmdsqsinvcut = Scalar(1.0) / rmdsqscut;
                Scalar rmdsqminvcut = Scalar(1.0) / rmdsqmcut;

                Scalar engt1cut = radprod * rmdsqsinvcut * A / Scalar(3.0);
                Scalar engt2cut = radprod * rmdsqminvcut * A / Scalar(3.0);
                Scalar engt3cut = slow::log(rmdsqscut * rmdsqminvcut) * A / Scalar(6.0);
                Scalar exp_valcut = fast::exp(-kappa * rmdscut);
                Scalar forcerepcut_divr = kappa * radprod * radsuminv * Z * exp_valcut / rcutt;
                pair_eng -= rcutt * forcerepcut_divr / kappa - engt1cut - engt2cut - engt3cut;
                }
            */
            //std::cout << "--new electrostatics--" << std::endl;
            //std::cout << "Z: " << 2RZ_over_kT << std::endl;
            //std::cout << "kappa: " << kappa << std::endl;
            //std::cout << "radsum: " << radsum << std::endl;
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
        return std::string("dlvoalt");
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
    unsigned int typei;       //!<~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;       //!<~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq;   //!< Stored rcutsq from the constructor
    Scalar cs_mM;    //!< kappa parameter extracted from the params passed to the constructor
    Scalar hb_nm;        //!< Z parameter extracted from the params passed to the constructor
    Scalar Zeta_mV;        //!< A parameter extracted from the params passed to the constructor
    Scalar radsum;   //!< radsum parameter extracted from the call to setDiameter
    Scalar radsub;   //!< radsub parameter extracted from the call to setDiameter
    Scalar radprod;  //!< radprod parameter extracted from the call to setDiameter
    Scalar radsumsq; //!< radsumsq parameter extracted from the call to setDiameter
    Scalar radsubsq; //!< radsubsq parameter extracted from the call to setDiameter
    Scalar delta;    //!< Diameter sum minus one
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_DLVO_ALT_H__
