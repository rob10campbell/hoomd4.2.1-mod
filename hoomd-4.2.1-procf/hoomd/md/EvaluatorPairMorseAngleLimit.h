// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __PAIR_EVALUATOR_MORSE_ANGLE_LIMIT_H__
#define __PAIR_EVALUATOR_MORSE_ANGLE_LIMIT_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h" // add vectors for optional position [RHEOINF]
#include "BondMap.h"
#include "hoomd/ParticleData.h"
#include <cmath> // for cosine calculation
#include <iostream>

/*! \file EvaluatorPairMorseAngleLimit.h
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

//extern BondDataManager manager;  // Declaration of the global BondDataManager object [RHEOINF]

namespace hoomd
    {
namespace md
    {

    //! Define functions that are used to calculate the angular repulsion potential 

    // create the Heaviside step function (2-1)
    // (this is for the angular repulsion potential)
    Scalar heaviside(Scalar r) {
        return (2 - r) < 0 ? 0.0 : 1.0;
    }

    // create the derivative of H(2-r), which is
    // a Dirac delta function d(H(2-r)/dr=-δ(r - 2)
    // (this is for the angular repulsion force)
    Scalar neg_dirac_delta(Scalar r, Scalar epsilon = 1e-5) {
        return std::abs(r - 2) < epsilon ? -1.0 / epsilon : 0.0;
    }

    // create the radial modulation function (Lambda function)
    // Λ(x) = x^-10 * [1-(0.5x)^10]^2 * heaviside(2-x)
    // (this is for the angular repulsion potential)
    Scalar lambda(Scalar x) {
        if (x == 0) return 0;  // Handle x == 0 separately to avoid division by zero

        Scalar h_lambda = heaviside(x);
        Scalar term = std::pow(0.5 * x, 10);
        return std::pow(x, -10) * std::pow(1 - term, 2) * h_lambda;
    }
            
    // and create it's derivative
    // Λ'(x) = -10x^-11 * [1-(0.5x)^10]^2 * heaviside(2-x)
    //            + r^-10 * ( 2[1-(r/2)^10]*(-5(r/2)^9 ) * heaviside(2-9)
    //               + r^-10 * [1-(r/2)^10]^2 * neg_dirac_delta
    // (this is for the angular repulsion force)
    Scalar lambda_prime(Scalar x) {
        if (x == 0) return 0;  // Handle x == 0 separately to avoid division by zero

        Scalar h_lambdaprime = heaviside(x);
        Scalar dirac_lambdaprime = neg_dirac_delta(x);

        Scalar halfterm_ten = std::pow(0.5 * x, 10);
        Scalar halfterm_nine = std::pow(0.5 * x, 9);

        Scalar second_term = std::pow(1 - halfterm_ten, 2);
        Scalar second_term_prime = -(5 * halfterm_nine);

        Scalar part_one = -10 * std::pow(x, -11) * second_term * h_lambdaprime;
        Scalar part_two = std::pow(x, -10) * (2 * second_term * second_term_prime) * h_lambdaprime;
        Scalar part_three = std::pow(x, -10) * second_term * dirac_lambdaprime;

        return part_one + part_two + part_three;
    }


//! Class for evaluating the Morse pair potential with angular repulsion in certain many-body angle arrangements
/*! <b>General Overview</b>

    See EvaluatorPairLJ.

    <b>MorseAngleLimit specifics</b>

    EvaluatorPairMorseAngleLimit evaluates the function:
    \f[ V_{\mathrm{Morse}}(r) = D_0 \left[ \exp \left(-2\alpha \left(r - r_0\right) \right)
                                           -2\exp \left(-\alpha \left(r-r_0\right) \right)  \right]
   \f]
*/
class EvaluatorPairMorseAngleLimit
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        Scalar r0;
        Scalar f_contact; //~ add f_contact param [RHEOINF]
        Scalar theta_bar;
        Scalar B;
        Scalar w;
        bool scaled_D0; //~ add scaled_D0 param [RHEOINF]

        //std::shared_ptr<BondMap> LTIME; // add pointer to Lifetime.h data

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        // CUDA memory hints
        void set_memory_hints() const { }
#endif

#ifndef __HIPCC__
        param_type() : D0(0), alpha(0), r0(0), f_contact(0), theta_bar(0), B(0), w(0), scaled_D0(false) { } //~ add f_contact and scaled_D0 params [RHEOINF]

        param_type(pybind11::dict v, bool managed = false)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            r0 = v["r0"].cast<Scalar>();
            f_contact = v["f_contact"].cast<Scalar>(); //~ add f_contact param [RHEOINF]
            theta_bar = v["theta_bar"].cast<Scalar>();
            B = v["B"].cast<Scalar>();
            w = v["w"].cast<Scalar>();
            this->scaled_D0 = scaled_D0; //~ add scaled_D0 param [RHEOINF]
            }

        param_type(Scalar d, Scalar a, Scalar r, Scalar f, Scalar tb, Scalar b, Scalar ww, bool sD, bool managed = false) //~ add f_contact and scaled_D0 params [RHEOINF]
            {
            D0 = d;
            alpha = a;
            r0 = r;
            f_contact = f; //~ add f_contact param [RHEOINF]
            theta_bar = tb;
            B = b;
            w = ww;
            scaled_D0 = sD; //~ add scaled_D0 param [RHEOINF]
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            v["r0"] = r0;
            v["f_contact"] = f_contact; //~ add f_contact param [RHEOINF]
            v["theta_bar"] = theta_bar;
            v["B"] = B;
            v["w"] = w;
            v["scaled_D0"] = scaled_D0; //~ add scaled_D0 param [RHEOINF] 
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
    DEVICE EvaluatorPairMorseAngleLimit(Scalar _rsq, Scalar _radcontact, unsigned int _pair_typeids[2], Scalar _rcutsq, const param_type& _params) //~add radcontact, pair_typeIDs [RHEOINF]
//        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), tag_i(0), tag_j(0), diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), f_contact(_params.f_contact), scaled_D0(_params.scaled_D0), bondMap() //~ add radcontact, diameters, f_contact, and scaled_D0 [RHEOINF]
        //: rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), tag_i(0), tag_j(0), diameter_i(0), diameter_j(0), p_i(0.0,0.0,0.0), p_j(0.0,0.0,0.0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), f_contact(_params.f_contact), theta_bar(_params.theta_bar), B(_params.B), w(_params.w), scaled_D0(_params.scaled_D0) //~ add radcontact, diameters, f_contact, and scaled_D0 [RHEOINF]
        : rsq(_rsq), radcontact(_radcontact), rcutsq(_rcutsq), tag_i(0), tag_j(0), diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), f_contact(_params.f_contact), theta_bar(_params.theta_bar), B(_params.B), w(_params.w), scaled_D0(_params.scaled_D0) //~ add radcontact, diameters, f_contact, and scaled_D0 [RHEOINF]
        {
        typei = _pair_typeids[0]; //~ add typei [RHEOINF]
        typej = _pair_typeids[1]; //~ add typej [RHEOINF] 
        }

    //! Whether the pair potential needs particle tags.
    HOSTDEVICE static bool needsTags()
        {
        return true;
        }
    //! Accept the optional tags
    /*! \param tag_i Tag of particle i
        \param tag_j Tag of particle j
    */
    HOSTDEVICE void setTags(unsigned int tagi, unsigned int tagj)
        {
        tag_i = tagi;
        tag_j = tagj;
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

    //!~ Whether the pair potential uses timestep. [RHEOINF]
    HOSTDEVICE static bool needsTimestep()
        {
        return true;
        }
    //! Accept the optional timestep
    /*! \param timestep the current timestep
    */
    HOSTDEVICE void setTimestep(uint64_t timestep)
        {
        currentTimestep = timestep;
        }

    //~ add i and j positions [RHEOINF] 
    DEVICE static bool needsIJPos()
        {
        return true;
        }
    //! Accept the optional position values
    /*! \param di position of particle i
        \param dj position of particle j
    */
    DEVICE void setIJPos(Scalar3 pi, Scalar3 pj)
        {
        p_i = pi;
        p_j = pj;
        }
    //~

    //!~ Whether the potential pair needs BoxDim info [RHEOINF]
    HOSTDEVICE static bool needsBox()
        {
        return true;
        }
    //! Accept the optional BoxDim structure
    /*! \param box the current box
    */
    HOSTDEVICE void setBox(const BoxDim box)
        {
        systembox = box;
        }


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

            //~ add contact force [RHEOINF]
            //~ check if contact force is provided [RHEOINF]
            if (f_contact != 0.0)
            {
                //~ if particles overlap (r < radsum) apply contact force
                if(r < radsum)force_divr = f_contact * (Scalar(1.0) - (r-radsum)) * pow((Scalar(0.50)*radsum),3) / r;

                else{
                    //~ calculate force as normal
                    force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;

                    //~ but still include contact force within 0.001 dist of colloid-colloid contact 
                    Scalar Del_max = Scalar(0.001); //~ 0.001 or 0.01
                    if(r<(radsum+Del_max))force_divr += f_contact * pow((Scalar(1.0) - (r-radsum)/Del_max), 3) * pow((Scalar(0.50)*radsum),3) / r;
                    } 
            
            }

            else {
                //~ calculate force as normal
                force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r;
                 }
            //~

            pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
            //~ force_divr = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r; //~ move this into overlap check [RHEOINF]


            //~ check BondMap [RHEOINF]
            //int currentData = manager.getData(); // if the data in BondMap is an int
            //std::cout << "CurrentData: " << currentData << std::endl;
            // ints from map
            //int value1 = manager.getData("key1");
            //int value2 = manager.getData("key2");
            //std::cout << "CurrentData: " << value1 << "," << value2 << std::endl;
            // bool and vaues from map
            BondDataManager& manager = BondDataManager::getInstance();
            //int value1;
            //int value2;
            //if (manager.getData(tag_i, value1)) {
            //    std::cout << "Value for key1: " << value1 << std::endl;
            //} else {
            //    std::cerr << "Key for tag " << tag_i << " not found" << std::endl;
            //}

            //if (manager.getData(tag_i, value1)) {
            //    int result = value1;
            //    std::cout << "tag_i (" << tag_i << "): " << result << std::endl;
            //}
            BondTrackData data;
            //const BoxDim box = this->m_pdata->getBox(); // get box for scaling
            Scalar3 rij_vector = p_j-p_i; // calculate r as a vector
            rij_vector = systembox.minImage(rij_vector); // scale vector within box

            // if bond map data exists
            if (manager.getData(tag_i, tag_j, data)) {
                // and if there are more than 1 bonds in the data set
                int num_bonds = manager.countTagIKeys(tag_i);
                if (num_bonds > 1) {
                    //std::cout << num_bonds << " bonds found, calculating many-body angular repulsion" << std::endl;
                    //std::cout << "tag_i, tag_j (" << tag_i << "," << tag_j<< ") " << num_bonds << " bonds" << std::endl;

                    // calculate the angular repulsion from each other existing bond

                    // find all data that includes tag_i
                    //std::vector<BondTrackData> bondDataVec = manager.getDataForTag(tag_i);
                    auto bondDataVec = manager.getDataForTag(tag_i);
                    //for (const auto& data : bondDataVec) {
                    for (const auto& pair : bondDataVec) {
                        const auto& key = pair.first;
                        const auto& data = pair.second;
  
                        // bond i-j does not interact with itself 
                        if (!((key.first == tag_i && key.second == tag_j) || (key.first == tag_j && key.second == tag_i))) { 
                            //for all other i-k pairs
                            Scalar3 pi = make_scalar3(0.0,0.0,0.0);
                            Scalar3 pk = make_scalar3(0.0,0.0,0.0);
                            if (key.first == tag_i) {
                                // r_prime is the vector from key.first to key.second
                                pi = make_scalar3(data.pos1_x, data.pos1_y, data.pos1_z);
                                pk = make_scalar3(data.pos2_x, data.pos2_y, data.pos2_z);
                            } else if (key.second == tag_i ){
                                // r_prime is the vector from key.second to key.first
                                pi = make_scalar3(data.pos2_x, data.pos2_y, data.pos2_z);
                                pk = make_scalar3(data.pos1_x, data.pos1_y, data.pos1_z);
                            }                  
                            Scalar3 r_prime = pk-pi; // calculate r_prime
                            r_prime = systembox.minImage(r_prime);
                            Scalar r_prime_mag = std::sqrt(r_prime.x * r_prime.x + r_prime.y * r_prime.y + r_prime.z * r_prime.z);

                            // calculate the vector dot product
                            Scalar dot_product = rij_vector.x * r_prime.x + rij_vector.y * r_prime.y + rij_vector.z * r_prime.z;

                            // calculate the angle in radians
                            Scalar cos_theta = dot_product / (r * r_prime_mag);

                            // make sure correct -1 <= cos_theta <= 1
                            if (cos_theta > 1.0)
                                cos_theta = 1.0;
                            if (cos_theta < -1.0)
                                cos_theta= -1.0;  

                            if (cos_theta != std::cos(theta_bar)) {
                                //std::cout << "Angular repulsion will be applied" << std::endl;

                                // calculate energy
                                Scalar lambda_r = lambda(r);
                                Scalar lambda_r_prime = lambda(r_prime_mag);
                                Scalar wsq_inv = (1 / std::pow(w, 2));
                                Scalar angular_repulsion_exp = fast::exp(-std::pow(cos_theta - std::cos(theta_bar), 2) * wsq_inv ); 
                                pair_eng += B * lambda_r * lambda_r_prime * angular_repulsion_exp;

                                // calculate force
                                Scalar lambdaprime_r = lambda_prime(r);
                                Scalar ddr_cos_theta = ( ((r_prime_mag * cos_theta) * (r * r_prime_mag)) - (dot_product * r_prime_mag) ) / std::pow(r * r_prime_mag , 2);
                                Scalar ddr_exp_contents = -2 * (cos_theta - std::cos(theta_bar)) * wsq_inv * ddr_cos_theta;
                                Scalar ddr_full_exp_term = angular_repulsion_exp * ddr_exp_contents;
                                Scalar first_force_term = lambdaprime_r * angular_repulsion_exp;
                                Scalar second_force_term = lambda_r * ddr_full_exp_term;
                                force_divr += B * lambda_r_prime * (first_force_term + second_force_term) / r;
                            }

                        }

                        //std::cout << "Key: (" << key.first << ", " << key.second << ")" << std::endl;
                        //std::cout << "formedTime=" << data.formedTime << "brokeTime=" << data.brokeTime << ", pos1=(" << data.pos1_x << ", "
                        //      << data.pos1_y << ", " << data.pos1_z << "), pos2=(" << data.pos2_x << ", "
                        //      << data.pos2_y << ", " << data.pos2_z << "), d1=" << data.d1 << ", d2=" << data.d2 << ", type1=" << data.type1 
                        //      << ", type2=" << data.type2 << std::endl;
                    }

                }
            //    std::cout << "Data found: formedTime=" << data.formedTime << ", brokeTime=" << data.brokeTime << ", pos1=(" << data.pos1_x << ", " << data.pos1_y << ", " << data.pos1_z << "), pos2=(" << data.pos2_x << ", " << data.pos2_y << ", " << data.pos2_z << "), d1=" << data.d1 << ", d2=" << data.d2 << ", type1" << data.type1 << ", type2=" << data.type2 << std::endl;
            //} else {
            //if (!manager.getData(tag_i, tag_j, data)) {
            //    std::cout << "Data not found for tags (" << tag_i << ", " << tag_j << ")" << std::endl;
            }
 
            //std::cout << "-----" << std::endl;
            //~



            if (energy_shift)
                {
                Scalar rcut = fast::sqrt(rcutsq);
                Scalar Exp_factor_cut = fast::exp(-alpha * (rcut - radsum));
                //Scalar Exp_factor_cut = fast::exp(-alpha * (rcut - r0));
                pair_eng -= D0 * Exp_factor_cut * (Exp_factor_cut - Scalar(2.0));
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
        return std::string("morse");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    //BondDataManager manager; // Implementation of BondDataManager object [RHEOINF]
    protected:
    const std::shared_ptr<ParticleData> m_pdata;
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar radcontact;//!< Stored contact-distance from the constructor [RHEOINF]
    unsigned int pair_typeids;//!< Stored pair typeIDs from the constructor [RHEOINF]
    unsigned int typei;//!<~ Stored typeID of particle i from the constructor [RHEOINF]
    unsigned int typej;//!<~ Stored typeID of particle j from the constructor [RHEOINF]
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    unsigned int tag_i;
    unsigned int tag_j;
    Scalar diameter_i;//!<~ add diameter_i [RHEOINF]
    Scalar diameter_j;//!<~ add diameter_j [RHEOINF]
    Scalar3 p_i;   //!< position of particle i
    Scalar3 p_j;    //!< position of particle i
    Scalar D0;     //!< Depth of the Morse potential at its minimum
    Scalar alpha;  //!< Controls width of the potential well
    Scalar r0;     //!< Offset, i.e., position of the potential minimum
    Scalar f_contact; //!< Contact force magnitude, for resolving overlap [RHEOINF]
    Scalar theta_bar; //!< allowed reference angle for multi-body bonds
    Scalar B;         //!< magnitude of the angular repulsion potential
    Scalar w;         //!< width the angular repulsion potential (about theta_bar)
    bool scaled_D0;   //!<~ on/off bool for scaling D0 by particle size [RHEOINF]
    //BondMap bondMap; //!< Reference to the bond map [RHEOINF]
    uint64_t currentTimestep;
    BoxDim systembox; //!< system box information passed from PotentialPair
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_ANGLE_LIMIT_H__
