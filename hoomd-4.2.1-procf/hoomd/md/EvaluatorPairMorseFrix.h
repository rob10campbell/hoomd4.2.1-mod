// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// $Id$
// $URL$

#ifndef __PAIR_EVALUATOR_MORSE_FRIX_H__
#define __PAIR_EVALUATOR_MORSE_FRIX_H__

#ifndef __HIPCC__
#include <string>
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h" // add vectors for optional position [RHEOINF]
#endif

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif
#include "hoomd/VectorMath.h"
#include <iostream>
/*! \file EvaluatorPairMorseFrix.h
    \brief Defines the Morse potential with frictional restrictions on rotation
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
class EvaluatorPairMorseFrix
    {
    public:
    struct param_type
        {
        Scalar D0;
        Scalar alpha;
        //Scalar r0;
        Scalar kn;
        Scalar eta_t;
        Scalar mu_t;

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

        //HOSTDEVICE param_type() : D0(0), alpha(0), r0(0), kn(0), eta_t(0), mu_t(0) { }
        HOSTDEVICE param_type() : D0(0), alpha(0), kn(0), eta_t(0), mu_t(0) { }
        //HOSTDEVICE param_type() : kn(0), eta_t(0), mu_t(0) { }

#ifndef __HIPCC__

        param_type(pybind11::dict v, bool managed)
            {
            D0 = v["D0"].cast<Scalar>();
            alpha = v["alpha"].cast<Scalar>();
            //r0 = v["r0"].cast<Scalar>();
            kn = v["kn"].cast<Scalar>();
            eta_t = v["eta_t"].cast<Scalar>();
            mu_t = v["mu_t"].cast<Scalar>();
            }

        pybind11::object toPython()
            {
            pybind11::dict v;
            v["D0"] = D0;
            v["alpha"] = alpha;
            //v["r0"] = r0;
            v["kn"] = kn;
            v["eta_t"] = eta_t;
            v["mu_t"] = mu_t;
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
        \param _kn
        \param _eta_t 
        \param _mu_t 
        \param _params Per type pair parameters of this potential
    */
    HOSTDEVICE EvaluatorPairMorseFrix(Scalar3& _dr,
                                     Scalar4& _quat_i,
                                     Scalar4& _quat_j,
                                     Scalar _rcutsq,
                                     const param_type& _params
                                     )
        : dr(_dr), rcutsq(_rcutsq), tag_i(0), tag_j(0), quat_i(_quat_i), quat_j(_quat_j),
          //mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), kn(_params.kn), eta_t(_params.eta_t), mu_t(_params.mu_t)
          mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), kn(_params.kn), eta_t(_params.eta_t), mu_t(_params.mu_t)
          //mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), kn(_params.kn), eta_t(_params.eta_t), mu_t(_params.mu_t)
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

    //~ Whether the pair potential uses velocity [RHEOINF]
    HOSTDEVICE static bool needsVel()
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

    // Accept optional velocity [RHEOINF]
    /* \param  v_tan velocity
       from AnisoPotentialPair.h
    */
    HOSTDEVICE void setVel(vec3<Scalar> U_t)
    {
    v_tan = U_t;
    }

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
        Scalar r = fast::sqrt(rsq);
        Scalar radsum = 0.5 * (diameter_i + diameter_j); //~ Add radsum from passed diameters
        Scalar h_ij = r - radsum;

        Scalar rcut = fast::sqrt(rcutsq);

        //std::cout << "Frix..." << std::endl;

        // if the particles do not interact
        if (rsq > rcutsq)
           {
           return false;
           }
 
        // Compute Morse 
        // /*         
        //Scalar r = fast::sqrt(rsq); 
        Scalar Exp_factor = fast::exp(-alpha * (r - radsum));

        //~ add contact force [RHEOINF]
        //~ check if contact force is provided [RHEOINF]
        Scalar f_contact = 0.0;
        Scalar3 Morse_force = make_scalar3(0, 0, 0);
        if (f_contact != 0.0)
        {
            //~ if particles overlap (r < radsum) apply contact force
            if(r < radsum)Morse_force = vec_to_scalar3(f_contact * (Scalar(1.0) - (r-radsum)) * pow((Scalar(0.50)*radsum),3) * rvec / r);

            else{
                //~ calculate force as normal
                Morse_force += vec_to_scalar3(Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) * rvec / r);

                //~ but still include contact force within 0.001 dist of colloid-colloid contact 
                Scalar Del_max = Scalar(0.001); //~ 0.001 or 0.01
                if(r<(radsum+Del_max))Morse_force += vec_to_scalar3(f_contact * pow((Scalar(1.0) - (r-radsum)/Del_max), 3) * pow((Scalar(0.50)*radsum),3) * rvec / r);
                } 
           
        }

        else {
            //~ calculate force as normal
            Morse_force += vec_to_scalar3(Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) * rvec / r);
             }
        //~
        pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
        //~ Morse_force = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) / r; //~ move this into overlap check [RHEOINF]
        //std::cout << D0 << "," << Exp_factor << "," << pair_eng << std::endl;            

        //Scalar Exp_factor = fast::exp(-alpha * (r - radsum));
        //f_scalar = Scalar(2.0) * D0 * alpha * Exp_factor * (Exp_factor - Scalar(1.0)) * rinv;
        //pair_eng = D0 * Exp_factor * (Exp_factor - Scalar(2.0));
        //*/

        Scalar range = 3/alpha;
        Scalar frix_cut = 0.5*range;
        Scalar del_min = 0.001;
        if (h_ij > del_min and h_ij < rcut) {
        // Computer Friction
        Scalar overlap = (1.0 - (h_ij / rcut)) * 0.1; // 0 at h_ij=r_cut, 0.1 at colloid-colloid contact (prevents force blow up)
        //Scalar overlap = 1 / (h_ij / rcut); // at rcut, overlap = 1, as particles get closer, the friction increases
        //vec3<Scalar> F_normal = -kn * overlap * rvec / r; //original: but - = repulsive ; -(j->i) AKA i->j
        //vec3<Scalar> F_normal = kn * overlap * rvec / r; //suggested: but + = attractive ; j->i
        vec3<Scalar> F_normal = kn * overlap * rvec / r; //suggested: but + = attractive ; j->i
        //vec3<Scalar> F_normal = Morse_force + (kn * overlap * rvec / r); //suggested: but + = attractive ; j->i
        //vec3<Scalar> Ft_damp = -eta_t * v_tan;
        vec3<Scalar> Ft_damp = eta_t * v_tan;
        //Scalar F_normal_mag = fast::sqrt(dot(F_normal,F_normal));
        Scalar F_normal_mag = fast::sqrt(dot(F_normal,F_normal)) + fast::sqrt(dot(Morse_force,Morse_force));
        Scalar Ft_damp_mag =  fast::sqrt(dot(Ft_damp,Ft_damp));
        vec3<Scalar> t = normalize(v_tan);
        //vec3<Scalar> friction_F = - std::min(mu_t*F_normal_mag,Ft_damp_mag) * t;
        //vec3<Scalar> friction_F = std::min(mu_t*F_normal_mag,Ft_damp_mag) * t;
        vec3<Scalar> friction_F = mu_t*F_normal_mag * t;
        //vec3<Scalar> friction_F = mu_t*F_normal_mag * t;
        force+= vec_to_scalar3(friction_F + F_normal);

        torque_i -= vec_to_scalar3(cross((diameter_i+0.5*overlap)*rvec/r,friction_F));
        torque_j -= vec_to_scalar3(cross((diameter_j+0.5*overlap)*rvec/r,friction_F));
        //torque_i += vec_to_scalar3(cross((diameter_i+0.5*overlap)*rvec/r,friction_F));
        //torque_j += vec_to_scalar3(cross((diameter_j+0.5*overlap)*rvec/r,friction_F));

        } // close correct r_cut = 3/kappa range

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
        return "morsefrix";
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
    Scalar D0;     //!< Depth of the Morse potential at its minimum
    Scalar alpha;  //!< Controls width of the potential well
    //Scalar r0;     //!< Offset, i.e., position of the potential minimum
    Scalar kn;
    Scalar eta_t;
    Scalar mu_t;
    vec3<Scalar> v_tan;
    // const param_type &params;   //!< The pair potential parameters
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_FRIX_H__
