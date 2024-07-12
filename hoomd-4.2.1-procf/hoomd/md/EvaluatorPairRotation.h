// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// $Id$
// $URL$

#ifndef __PAIR_EVALUATOR_ROTATION_H__
#define __PAIR_EVALUATOR_ROTATION_H__

#ifndef __HIPCC__
#include <string>
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h" // add vectors for optional position [RHEOINF]
#include "hoomd/md/RotationMap.h"
#endif

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif
#include "hoomd/VectorMath.h"
#include <iostream>
/*! \file EvaluatorPairRotation.h
    \brief Defines the angular rotation potential
*/

// need to declare these class methods with __device__ qualifiers when building
// in nvcc.  HOSTDEVICE is __host__ __device__ when included in nvcc and blank
// when included into the host compiler
#ifdef __HIPCC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif

namespace hoomd
    {
namespace md
    {
class EvaluatorPairRotation
    {
    public:
    struct param_type
        {
        //Scalar D0;
        //Scalar alpha;
        //Scalar r0;
        Scalar K; //! The energy oscillation magnitude encountered during rotation 
        Scalar n; //! The energy oscillation frequency (# of peaks of height 2K)

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

        //! Load dynamic data members into shared memory and increase pointer
        /*! \param ptr Pointer to load data to (will be incremented)
            \param available_bytes Size of remaining shared memory
            allocation
        */
        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

        //HOSTDEVICE param_type() : D0(0), alpha(0), r0(0), K(0), n(0) { }
        HOSTDEVICE param_type() : K(0), n(0) { }

#ifndef __HIPCC__

        param_type(pybind11::dict v, bool managed)
            {
            //D0 = v["D0"].cast<Scalar>();
            //alpha = v["alpha"].cast<Scalar>();
            //r0 = v["r0"].cast<Scalar>();
            K = v["K"].cast<Scalar>();
            n = v["n"].cast<Scalar>();
            }

        pybind11::object toPython()
            {
            pybind11::dict v;
            //v["D0"] = D0;
            //v["alpha"] = alpha;
            //v["r0"] = r0;
            v["K"] = K;
            v["n"] = n;
            return std::move(v);
            }

#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    struct shape_type
        {
        vec3<Scalar> mu;
        //! Load dynamic data members into shared memory and increase pointer
        /*! \param ptr Pointer to load data to (will be incremented)
            \param available_bytes Size of remaining shared memory allocation
        */
        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

        HOSTDEVICE shape_type() : mu {0, 0, 0} { }

#ifndef __HIPCC__

        shape_type(vec3<Scalar> mu_, bool managed = false) : mu(mu_) { }

        shape_type(pybind11::object mu_obj, bool managed)
            {
            auto mu_ = (pybind11::tuple)mu_obj;
            mu = vec3<Scalar>(mu_[0].cast<Scalar>(), mu_[1].cast<Scalar>(), mu_[2].cast<Scalar>());
            }

        pybind11::object toPython()
            {
            return pybind11::make_tuple(mu.x, mu.y, mu.z);
            }
#endif // __HIPCC__

#ifdef ENABLE_HIP
        //! Attach managed memory to CUDA stream
        void set_memory_hint() const { }
#endif
        };

    //! Constructs the pair potential evaluator
    /*! \param _dr Displacement vector between particle centers of mass
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _quat_i Quaternion of i^{th} particle
        \param _quat_j Quaternion of j^{th} particle
        \param _K energy oscillation magnitude
        \param _n energy oscillation frequency
        \param _params Per type pair parameters of this potential
    */
    HOSTDEVICE EvaluatorPairRotation(Scalar3& _dr,
                                     Scalar4& _quat_i,
                                     Scalar4& _quat_j,
                                     Scalar _rcutsq,
                                     const param_type& _params
                                     )
        : dr(_dr), rcutsq(_rcutsq), tag_i(0), tag_j(0), quat_i(_quat_i), quat_j(_quat_j),
          //mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), K(_params.K), n(_params.n), rotationMap()
          mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), K(_params.K), n(_params.n), rotationMap()
        {
        }

    //! Whether the pair potential uses shape.
    HOSTDEVICE static bool needsShape()
        {
        return true;
        }

    //! Whether the pair potential needs particle tags.
    HOSTDEVICE static bool needsTags()
        {
        return true;
        }
        
    //! don't need diameter
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

    //! whether pair potential requires charges
    HOSTDEVICE static bool needsCharge()
        {
        return false;
        }

    //!~ Whether the pair potential uses typeid. [RHEOINF]
    HOSTDEVICE static bool needsTypes()
        {
        return true;
        }
    //~

    //!~ Whether the pair potential uses timestep. [RHEOINF]
    HOSTDEVICE static bool needsTimestep()
        {
        return true;
        }
    //~


    /// Whether the potential implements the energy_shift parameter
    HOSTDEVICE static bool constexpr implementsEnergyShift()
        {
        return false;
        }


    /*! \param shape_i Shape of particle i
        \param shape_j Shape of particle j
    */
    HOSTDEVICE void setShape(const shape_type* shapei, const shape_type* shapej)
        {
        mu_i = shapei->mu;
        mu_j = shapej->mu;
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

    //! Accept the optional charge values
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    HOSTDEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Accept the optional types
    /*! \param type_i typeID of particle i
        \param type_j typeID of particle j
    */
    HOSTDEVICE void setTypes(unsigned int typei, unsigned int typej) 
    {
    type_i = typei;
    type_j = typej;
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
    /*! \param force Output parameter to write the computed force.
        \param pair_eng Output parameter to write the computed pair energy.
        \param energy_shift If true, the potential must be shifted so that
            V(r) is continuous at the cutoff.
        \param torque_i The torque exerted on the i^th particle.
        \param torque_j The torque exerted on the j^th particle.
        \return True if they are evaluated or false if they are not because
            we are beyond the cutoff.
    */
    HOSTDEVICE bool evaluate(Scalar3& force,
                             Scalar& pair_eng,
                             bool energy_shift,
                             Scalar3& torque_i,
                             Scalar3& torque_j)
        {
        vec3<Scalar> rvec(dr);
        Scalar rsq = dot(rvec, rvec);

        std::cout << "Rotation Potential..." << std::endl;

        // if the particles do not interact
        if (rsq > rcutsq)
           {
           // check if the bond broke
           std::cout << "rsq > rcutsq" << std::endl;
           std::cout << "rsq = " << rsq << "rcutsq = " << rcutsq << std::endl;
           BondInfo bond_info;
           if (rotationMap.getBondInfo(tag_i, tag_j, bond_info))
           {
             // Bond existed before this, mark it as broken
             rotationMap.addBondBroke(currentTimestep, tag_i, tag_j);
           }

           return false;
           }

        //~ Add radsum from passed diameters [RHEOINF] 
        //Scalar radsum = 0.5 * (diameter_i + diameter_j);
 
        // Compute the current orientation vectors in the space frame
        // Assuming the initial orientation vector is (1, 0, 0) (x-axis)
        vec3<Scalar> initial_orientation_i(1, 0, 0);
        vec3<Scalar> initial_orientation_j(1, 0, 0);

        vec3<Scalar> e_i = rotate(quat<Scalar>(quat_i), mu_i);
        vec3<Scalar> e_j = rotate(quat<Scalar>(quat_j), mu_j);

        vec3<Scalar> f;
        vec3<Scalar> t_i;
        vec3<Scalar> t_j;
        Scalar e = Scalar(0.0);

        Scalar rinv = fast::rsqrt(rsq);
        Scalar r = Scalar(1.0) / rinv;
        //Scalar r2inv = Scalar(1.0) / rsq;

        // Compute Morse 
        //Scalar Exp_factor = fast::exp(-alpha * (r - radsum));
        //f_scalar = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) * rinv;
        //pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));


        // compute angles
        Scalar orientation_i = dot(e_i, rvec) / (sqrt(dot(e_i, e_i)) * r);
        Scalar orientation_j = dot(e_j, rvec) / (sqrt(dot(e_j, e_j)) * r);
        Scalar theta_i = acos(orientation_i);
        Scalar theta_j = acos(orientation_j);

        Scalar torsional_orientation_ij = dot(e_i, e_j) / (sqrt(dot(e_i, e_i)) * sqrt(dot(e_j, e_j)));
        Scalar gamma_ij = acos(torsional_orientation_ij);

        ////~ reference add bond_calc from PotentialPairDPDThermo [RHEOINF]
        //if(m_bond_calc)
        //  {
        //  if(typei && typej) // if both are NOT zero (solvents are type zero)
        //    {
        //    Scalar rsq_root = fast::sqrt(rsq) - Scalar(0.5)*(h_diameter.data[i]+h_diameter.data[j]);
        //    if(rsq_root < Scalar(0.10)) // assumes the cut-off is 0.1
        //      {
        //      if(rsq_root < Scalar(0.08))
        //        {
        //        this->LTIME->Bond_check[bond_index] = 2;
        //        }
        //      else
        //        {
        //        this->LTIME->Bond_check[bond_index] = 1;
        //        }

        std::cout << "Checking RotationMap" << std::endl;

        BondInfo bond_info;
        if (!rotationMap.getBondInfo(tag_i, tag_j, bond_info))
        {
            std::cout << "New bond found" << std::endl;
            std::cout << "Adding bond to the RotationMap" << std::endl;
            // Bond is new, add it to the map
            rotationMap.addBondFormed(currentTimestep, tag_i, tag_j, type_i, type_j, theta_i, theta_j, gamma_ij);
            rotationMap.getBondInfo(tag_i, tag_j, bond_info);
            std::cout << bond_info.timestep << std::endl;
            std::cout << bond_info.id1 << std::endl;
            std::cout << bond_info.id2 << std::endl;
            std::cout << bond_info.typei << std::endl;
            std::cout << bond_info.typej << std::endl;
            std::cout << bond_info.formedTimestep << std::endl;
            std::cout << bond_info.theta_i_0 << std::endl;
            std::cout << bond_info.theta_j_0 << std::endl;
            std::cout << bond_info.gamma_ij_0 << std::endl;
            std::cout << bond_info.broken << std::endl;
            std::cout << bond_info.brokeTimestep << std::endl;

        }
        else
            {

            std::cout << "Extracting RotationMap" << std::endl;
            rotationMap.getBondInfo(tag_i, tag_j, bond_info);
            std::cout << bond_info.id1 << std::endl;
            std::cout << bond_info.id2 << std::endl;
            std::cout << bond_info.typei << std::endl;
            std::cout << bond_info.typej << std::endl;
            std::cout << bond_info.formedTimestep << std::endl;
            std::cout << bond_info.theta_i_0 << std::endl;
            std::cout << bond_info.theta_j_0 << std::endl;
            std::cout << bond_info.gamma_ij_0 << std::endl;
            std::cout << bond_info.broken << std::endl;
            std::cout << bond_info.brokeTimestep << std::endl;
            // Bond exists, use the stored values
            Scalar theta_i_0 = bond_info.theta_i_0;
            Scalar theta_j_0 = bond_info.theta_j_0;
            Scalar gamma_ij_0 = bond_info.gamma_ij_0;

            std::cout << "Calculating amount of rotation since last step" << std::endl;
            Scalar del_theta_i = theta_i - theta_i_0;
            Scalar del_theta_j = theta_j - theta_j_0;
            Scalar del_gamma_ij = gamma_ij - gamma_ij_0;

            // Compute rotation potential
            Scalar U_rot = (K / r) * (3 - cos(n * del_theta_i) - cos(n * del_theta_j) - cos(n * del_gamma_ij));

            std::cout << "Pair U_rot: " << U_rot << std::endl;

            // Compute the force F_rot            
            //
            //f_scalar = 
            //vec3<Scalar> f = f_scalar * (dr / r);
            //

            e += U_rot;

            force = vec_to_scalar3(f);
            torque_i = vec_to_scalar3(t_i);
            torque_j = vec_to_scalar3(t_j);
            pair_eng = e;

            }

        std::cout << "Writing current Rotationmap to BondHistory file" << std::endl;
        rotationMap.writeBondHistory(currentTimestep, tag_i, tag_j);

        return true;
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
    //! Get the name of the potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return "rotation";
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar3 dr;             //!< Stored vector pointing between particle centers of mass
    Scalar rcutsq;          //!< Stored rcutsq from the constructor
    unsigned int tag_i;
    unsigned int tag_j;
    unsigned int type_i;
    unsigned int type_j;
    Scalar4 quat_i, quat_j; //!< Stored quaternion of ith and jth particle from constructor
    vec3<Scalar> mu_i;      /// Magnetic moment for ith particle
    vec3<Scalar> mu_j;      /// Magnetic moment for jth particle
    Scalar diameter_i;//!<~ add diameter_i [RHEOINF]
    Scalar diameter_j;//!<~ add diameter_j [RHEOINF]
    //Scalar D0;     //!< Depth of the Morse potential at its minimum
    //Scalar alpha;  //!< Controls width of the potential well
    //Scalar r0;     //!< Offset, i.e., position of the potential minimum
    Scalar K;
    Scalar n;
    vec3<Scalar> e_i, e_j;
    RotationMap rotationMap; //!< Reference to the rotation map
    uint64_t currentTimestep;
    // const param_type &params;   //!< The pair potential parameters
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_ROTATION_H__
