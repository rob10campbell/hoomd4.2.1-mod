// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

//

#ifndef __EVALUATOR_ANGULAR_REPULSION__
#define __EVALUATOR_ANGULAR_REPULSION__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorAngularRepulsion.h
    \brief Defines the evaluator class for the three-body AngularRepulsion potential
*/

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
//! Class for evaluating the AngularRepulsion three-body potential
class EvaluatorAngularRepulsion
    {
    public:
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        Scalar r0;
        //Scalar f_contact
        Scalar theta_bar;
        Scalar B;
        Scalar w;
        bool scaled_D0;

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), r0(0), theta_bar(0), B(0), w(0), scaled_D0(false) { }

        param_type(pybind11::dict v)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            r0 = v["r0"].cast<Scalar>();
            //f_contact = v["f_contact"].cast<Scalar>();
            theta_bar = v["theta_bar"].cast<Scalar>();
            B = v["B"].cast<Scalar>();
            w = v["w"].cast<Scalar>();
            this->scaled_D0 = scaled_D0; 
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["r0"] = r0;
            //v["f_contact"] = f_contact;
            v["theta_bar"] = theta_bar;
            v["B"] = B;
            v["w"] = w;
            v["scaled_D0"] = scaled_D0;
            return v;
            }
#endif
        } __attribute__((aligned(16)));

    //! Constructs the evaluator
    /*! \param _rij_sq Squared distance between particles i and j
        \param _rcutsq Squared distance at which the potential goes to zero
        \param _params Per type-pair parameters for this potential
    */
    DEVICE EvaluatorAngularRepulsion(Scalar _rij_sq,
                             Scalar _rcutsq,
                             Scalar _diam_i, Scalar _diam_j, Scalar _diam_k,
                             const param_type& _params) // here it receives also r cutoff
        : rij_sq(_rij_sq), rcutsq(_rcutsq), diam_i(_diam_i), diam_j(_diam_j), diam_k(_diam_k), D0(_params.D0), alpha(_params.alpha), r0(_params.r0),
          theta_bar(_params.theta_bar), B(_params.B), w(_params.w), scaled_D0(_params.scaled_D0)
        {
        }

    //! Set the square distance between particles i and j
    DEVICE void setRij(Scalar rsq)
        {
        rij_sq = rsq;
        }

    //! Set the square distance between particles i and k
    DEVICE void setRik(Scalar rsq)
        {
        rik_sq = rsq;
        }

    //! Check if the potential needs the diameters of i and j
    DEVICE static bool needsDiDj()
        {
        return true;
        }

    //! iet the diameters for particles i and j
    DEVICE void setDiDj(Scalar di, Scalar dj)
        {
        diameter_i = di;
        diameter_j = dj;
        }

    //! Check if the potential needs the rij vector
    DEVICE static bool needsRijVec()
        {
        return true;
        }

    //! Set the vector between particles i and j
    DEVICE void setRijVec(Scalar3 _rij_vec)
        {
        rij_vec = _rij_vec;
        }

    //! Check if the potential needs diameter k
    DEVICE static bool needsDk()
        {
        return true;
        }

    //! Set the diameter for particle k
    DEVICE void setDk(Scalar dk)
        {
        diameter_k = dk;
        }

    //! Check if the potential needs the rik vector
    DEVICE static bool needsRikVec()
        {
        return true;
        }

    //! Set the vector between particles i and j
    DEVICE void setRikVec(Scalar3 _rik_vec)
        {
        rik_vec = _rik_vec;
        }

    //! This is a pure pair potential
    DEVICE static bool hasPerParticleEnergy()
        {
        return false;
        }

    //! We do not need chi
    DEVICE static bool needsChi()
        {
        return false;
        }

    //! We have ik-forces
    DEVICE static bool hasIkForce()
        {
        return true;
        }

    //! The AngularRepulsion potential needs the bond angle
    DEVICE static bool needsAngle()
        {
        return true;
        }

    //! Set the bond angle value
    //! \param _cos_th Cosine of the angle between ij and ik
    DEVICE void setAngle(Scalar _cos_th)
        {
        cos_th = _cos_th;
        }

    //! Check whether a pair of particles experience angular repulsion
    DEVICE bool areInteractive()
        {
        if (needsDk() and needsDiDj())
            {
            Scalar radsum_ik = Scalar(0.5) * (diameter_i + diameter_k);
            Scalar h_ik = fast::sqrt(rik_sq) - radsum_ik;
            Scalar rcut_ik = fast::sqrt(rcutsq);
            return (h_ik < rcut_ik) && (B != Scalar(0.0));
            }
        else
            {
            return (rik_sq < rcutsq) && (B != Scalar(0.0));
            }
        }

    //! Evaluate the repulsive and attractive terms of the force
    //DEVICE void evalRepulsiveAndAttractive(Scalar& invratio, Scalar& invratio2) { }
    DEVICE bool evalRepulsiveAndAttractive(Scalar& invratio, Scalar& invratio2)
        {
        if (needsDiDj())
            {
            invratio = 0.0;
            invratio2 = 0.0;
            Scalar radsum_ij = Scalar(0.5) * (diameter_i + diameter_j);
            Scalar h_ij = fast::sqrt(rij_sq) - radsum_ij;
            Scalar rcut_ij = fast::sqrt(rcutsq);
            return (h_ij < rcut_ij); //&& (B != Scalar(0.0));
            }
        else
            {
            invratio = 0.0;
            invratio2 = 0.0;
            return (rij_sq < rcutsq); //&& (epsilon_dev != Scalar(0.0)))
            }
        }
    //    if ((rij_sq < rcutsq) //&& (epsilon_dev != Scalar(0.0)))
    //        {
    //        // compute rij
    //        Scalar rij = fast::sqrt(rij_sq);

    //        // compute the power of the ratio
    //        invratio = fast::pow(sigma_dev / rij, n_dev);
    //        invratio2 = invratio * invratio;

    //        return true;
    //        }
    //    else
    //        return false;
    //    }

    //! We do not have to evaluate chi
    DEVICE void evalChi(Scalar& chi) { }

    //! We don't have a scalar ij contribution
    DEVICE void evalPhi(Scalar& phi) { }

    //! Evaluate the force and potential energy due to ij interactions
    DEVICE void evalForceij(Scalar invratio, // not used
                            Scalar invratio2, // not used
                            Scalar chi,  // not used
                            Scalar phi,  // not used
                            Scalar& bij, // not used
                            Scalar& force_divr,
                            Scalar& potential_eng)
        {
        // compute the ij force
        // the force term includes rij_sq^-1 from the derivative over the distance and a factor 0.5
        // to compensate for double countings

        //Scalar r_ij = fast::sqrt(rij_sq)
        diameter_i = diam_i;
        diameter_j = diam_j;
        diameter_k = diam_k;

        //~ Add radsum from passed diameters [RHEOINF]  
        //Scalar radsum_ij = Scalar(0.0);
        //if (needsDiDj())
        //  {
          Scalar radsum_ij = Scalar(0.5) * (diameter_i + diameter_j); 
        //  }

        //Scalar h_ij = r_ij - radsum_ij

        //std::cout << "Forceij; radsum_ij: " << radsum_ij << std::endl;

        //~ Scale attraction strength by particle size is scaled_D0 is true
        if (scaled_D0)
          {
          D0 = D0 * (Scalar(0.5)*radsum_ij);
          }   
        //~ 

        // compute the force divided by r in force_divr
        if (rij_sq < rcutsq)
           {
           Scalar r = fast::sqrt(rij_sq); 
           Scalar Exp_factor = fast::exp(-alpha * (r - radsum_ij));
           //force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) * Scalar(2.0) / rij_sq;
           //force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;
           force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r; 

           // compute the potential energy
           potential_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
        //   return true;
           }
        //else
        //   return false;
        }


    DEVICE void evalSelfEnergy(Scalar& energy, Scalar phi) { }

    //! Evaluate the forces due to ijk interactions
    DEVICE bool evalForceik(Scalar ijinvratio, // not used
                            Scalar& angle_potential, //Scalar ijinvratio2, // replaced with energy
                            Scalar chi, // not used
                            Scalar bij, // not used
                            Scalar3& force_divr_ij,
                            Scalar3& force_divr_ik)
        {

        diameter_i = diam_i;
        diameter_j = diam_j;
        diameter_k = diam_k;

        Scalar3 rij_vector = make_scalar3(0.0,0.0,0.0);
        Scalar3 rik_vector = make_scalar3(0.0,0.0,0.0);

        //Scalar radsum_ij = Scalar(0.0);
        //Scalar radsum_ik = Scalar(0.0);
        //if (needsDk() and needsDiDj())
        //   {
           Scalar radsum_ij = Scalar(0.5) * (diameter_i + diameter_j); 
           Scalar radsum_ik = Scalar(0.5) * (diameter_i + diameter_k); 
        //   }


        //std::cout << "Forceik; radsum_ij: " << radsum_ij << std::endl;
        //std::cout << "Forceik; radsum_ik: " << radsum_ij << std::endl;

        Scalar h_ik = fast::sqrt(rik_sq) - radsum_ik;
        Scalar rcut_ik = fast::sqrt(rcutsq);

        if (h_ik < rcut_ik)
            {
            // For compatibility with Tersoff I get 3d vectors in, but I need only to calculate
            // their modulus and I store it in the x component
            //Scalar force_divr_ij = IN_force_divr_ij.x;
            //Scalar force_divr_ik = IN_force_divr_ik.x;
   
            // compute rij, rik, rcut
            Scalar rij = fast::sqrt(rij_sq);
            Scalar rik = fast::sqrt(rik_sq);

            // get sin(theta)
            Scalar sin_th = std::sqrt(1.0 - cos_th * cos_th);
            if (sin_th < Scalar(0.001))
               sin_th = Scalar(0.001);
            sin_th = 1.0 / sin_th;

            Scalar Lambda_ij = Scalar(0.0);
            Scalar Lambda_ik = Scalar(0.0);
            //if (needsDk() and needsDiDj())
            //    {
            //    Lambda_ij = pow(rij/radsum_ij, -10) * pow(1 - pow(0.5 * (rij / radsum_ij), 10), 2);
            //    Lambda_ik = pow(rik/radsum_ik, -10) * pow(1 - pow(0.5 * (rik / radsum_ik), 10), 2);
            //    }
            //else
            //    {
            //    Lambda_ij = pow(rij, -10) * pow(1 - pow(0.5 * (rij), 10), 2);
            //    Lambda_ik = pow(rik, -10) * pow(1 - pow(0.5 * (rik), 10), 2);
            //    }
            Scalar term_ij = std::pow(0.5 * rij, 10);
            Lambda_ij = std::pow(rij, -10) * std::pow(1 - term_ij, 2);

            Scalar halfterm_ten_ij = std::pow(0.5 * rij, 10);
            Scalar halfterm_nine_ij = std::pow(0.5 * rij, 9);
            Scalar second_term_ij = std::pow(1 - halfterm_ten_ij, 2);
            Scalar second_term_prime_ij = -(5 * halfterm_nine_ij);
            Scalar part_one_ij = -10 * std::pow(rij, -11) * second_term_ij; 
            Scalar part_two_ij = std::pow(rij, -10) * (2 * second_term_ij * second_term_prime_ij);
            Scalar Lambda_prime_ij = part_one_ij + part_two_ij;


            Scalar term_ik = std::pow(0.5 * rik, 10);
            Lambda_ik = std::pow(rik, -10) * std::pow(1 - term_ik, 2);

            Scalar halfterm_ten_ik = std::pow(0.5 * rik, 10);
            Scalar halfterm_nine_ik = std::pow(0.5 * rik, 9);
            Scalar second_term_ik = std::pow(1 - halfterm_ten_ik, 2);
            Scalar second_term_prime_ik = -(5 * halfterm_nine_ik);
            Scalar part_one_ik = -10 * std::pow(rik, -11) * second_term_ik; 
            Scalar part_two_ik = std::pow(rik, -10) * (2 * second_term_ik * second_term_prime_ik);
            Scalar Lambda_prime_ik = part_one_ik + part_two_ik;

            //Scalar exp_term = exp(-pow((cos_th - cos(theta_bar)) / w, 2));

            //Scalar vab = (-2 * B * sin_th / pow(w, 2)) * Lambda_ij * Lambda_ik * exp_term * (cos_th - cos(theta_bar));
            //Scalar a11 = vab * cos_th / (rij * rij);
            //Scalar a12 = -vab / (rij * rik);
            //Scalar a22 = vab * cos_th / (rik * rik);

            Scalar wsq_inv = std::pow(w, -2);
            Scalar angular_repulsion_exp = std::exp(-std::pow(cos_th - std::cos(theta_bar), 2) * wsq_inv );

            if (needsRijVec() and needsRikVec())
               {
               rij_vector = rij_vec;
               rik_vector = rik_vec;
               }

            Scalar dot_product =  dot_product = rij_vector.x * rik_vector.x + rij_vector.y * rik_vector.y + rij_vector.z * rik_vector.z;
            Scalar ddr_cos_th = ( ((rik * cos_th) * (rij * rik)) - (dot_product * rik) ) / std::pow(rij * rik , 2);
            Scalar ddr_exp_contents = -2 * (cos_th - std::cos(theta_bar)) * wsq_inv * ddr_cos_th;
            Scalar ddr_full_exp_term = angular_repulsion_exp * ddr_exp_contents;
            Scalar first_force_term = Lambda_prime_ij * angular_repulsion_exp;
            Scalar second_force_term = Lambda_ij * ddr_full_exp_term;
            Scalar angle_force_theta = B * Lambda_ik * (first_force_term + second_force_term);


            Scalar a11 = angle_force_theta * cos_th / (rij * rij); // force on a from i
            //Scalar a12_ria = - angle_force_theta / ria_mag; // force on i from a
            //Scalar a12_rib = - rib_angle_force / rib_mag; // force on i from b
            Scalar a12 = -angle_force_theta / (rij * rik);
            Scalar a22 = angle_force_theta * cos_th / (rik * rik); // force on b from i


            // assign the ij forces
            force_divr_ij.x = (a11 * rij_vector.x + a12 * rik_vector.x) / rij; 
            force_divr_ij.y = (a11 * rij_vector.y + a12 * rik_vector.y) / rij;
            force_divr_ij.z = (a11 * rij_vector.z + a12 * rik_vector.z) / rij;
            // assign the ik forces
            force_divr_ik.x = (a22 * rij_vector.x + a12 * rik_vector.x) / rik;
            force_divr_ik.y = (a22 * rij_vector.y + a12 * rik_vector.y) / rik;
            force_divr_ik.z = (a22 * rij_vector.z + a12 * rik_vector.z) / rik;


            //angle_potential = (B * Lambda_ij * Lambda_ik * exp_term)/3;

            angle_potential = B * Lambda_ij * Lambda_ik * angular_repulsion_exp;


            return true;
            }
        else
            return false;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("angularrepulsion");
        }
#endif

    static const bool flag_for_AngularRepulsion = true;
    static const bool flag_for_RevCross = false;

    protected:
    Scalar rij_sq; //!< Stored rij_sq from the constructor
    Scalar rik_sq; //!< Stored rik_sq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar diam_i;
    Scalar diam_j;
    Scalar diam_k;
    Scalar cos_th;      //!< Cosine of the angle between rij and rik
    Scalar diameter_i;  //!< Diameter of particle i
    Scalar diameter_j;  //!< Diameter of particle j
    Scalar3 rij_vec;    //!< xyz-vector for distance between i and j
    Scalar diameter_k;  //!< Diameter of particle k
    Scalar3 rik_vec;    //!< xyz-vector for distance between i and k
    Scalar D0;          //!< Depth of the Morse potential at its minimum
    Scalar alpha;       //!< Controls width of the potential well
    Scalar r0;          //!< Offset, i.e., position of the potential minimum
    //Scalar f_contact;   //!< Contact force magnitude, for resolving overlap [RHEOINF]
    Scalar theta_bar;   //!< allowed reference angle for multi-body bonds
    Scalar B;           //!< magnitude of the angular repulsion potential
    Scalar w;           //!< width the angular repulsion potential (about theta_bar)
    bool scaled_D0;     //!<~ on/off bool for scaling D0 by particle size [RHEOINF]
    };

    } // end namespace md
    } // end namespace hoomd

#endif
