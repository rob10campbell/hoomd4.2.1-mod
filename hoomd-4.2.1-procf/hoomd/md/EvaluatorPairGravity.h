// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// $Id$
// $URL$

#ifndef __PAIR_EVALUATOR_GRAVITY_H__
#define __PAIR_EVALUATOR_GRAVITY_H__

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
/*! \file EvaluatorPairGravity.h
    \brief Defines the gravitational potential in the Y-direction
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
class EvaluatorPairGravity
    {
    public:
    struct param_type
        {
        //Scalar D0;
        //Scalar alpha;
        //Scalar range; //[AJ 2024]
        //Scalar r0;
        //Scalar kn;
        //Scalar eta_n; //[AJ 2024]
        //Scalar eta_t;
        //Scalar mu_t;
        Scalar F_G;

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
        //HOSTDEVICE param_type() : D0(0), alpha(0), kn(0), eta_t(0), mu_t(0) { }
        HOSTDEVICE param_type() : F_G(0) { }
        //HOSTDEVICE param_type() : range(0), kn(0), eta_n(0), eta_t(0), mu_t(0) { }  //[AJ 2024]
        //HOSTDEVICE param_type() : kn(0), eta_t(0), mu_t(0) { }

#ifndef __HIPCC__

        param_type(pybind11::dict v, bool managed)
            {
            F_G = v["F_G"].cast<Scalar>();
            //D0 = v["D0"].cast<Scalar>();
            //alpha = v["alpha"].cast<Scalar>();
            //range = v["range"].cast<Scalar>();  //[AJ 2024]
            //r0 = v["r0"].cast<Scalar>();
            //kn = v["kn"].cast<Scalar>();
            //eta_n = v["eta_n"].cast<Scalar>();  //[AJ 2024]
            //eta_t = v["eta_t"].cast<Scalar>();
            //mu_t = v["mu_t"].cast<Scalar>();
            }

        pybind11::object toPython()
            {
            pybind11::dict v;
            v["F_G"] = F_G;
            //v["D0"] = D0;
            //v["alpha"] = alpha;
            //v["range"] = range;  //[AJ 2024]
            //v["r0"] = r0;
            //v["kn"] = kn;
            //v["eta_n"] = eta_n;  //[AJ 2024]
            //v["eta_t"] = eta_t;
            //v["mu_t"] = mu_t;
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
        \param _eta_n
        \param _eta_t 
        \param _mu_t 
        \param _params Per type pair parameters of this potential
    */
    HOSTDEVICE EvaluatorPairGravity(Scalar3& _dr,
                                     Scalar4& _quat_i,
                                     Scalar4& _quat_j,
                                     Scalar _rcutsq,
                                     const param_type& _params
                                     )
        : dr(_dr), rcutsq(_rcutsq), tag_i(0), tag_j(0), quat_i(_quat_i), quat_j(_quat_j),
          //mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), range(_params.range), kn(_params.kn), eta_n(_params.eta_n), eta_t(_params.eta_t), mu_t(_params.mu_t) // [AJ2024]
          //mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), D0(_params.D0), alpha(_params.alpha), r0(_params.r0), kn(_params.kn), eta_t(_params.eta_t), mu_t(_params.mu_t)
          mu_i {0, 0, 0}, mu_j {0, 0, 0}, diameter_i(0), diameter_j(0), 
          F_G(_params.F_G)
          //D0(_params.D0), alpha(_params.alpha), kn(_params.kn), eta_t(_params.eta_t), mu_t(_params.mu_t)
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
    //HOSTDEVICE static bool needsVeltan()
    HOSTDEVICE static bool needsVel()
       {
       return true;
       }
    //~

    // //~ Whether the pair potential uses velocity [AJ2024]
    // HOSTDEVICE static bool needsVelnorm()
    //    {
    //    return true;
    //    }
    // //~

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
    //HOSTDEVICE void setVeltan(vec3<Scalar> U_t)
    HOSTDEVICE void setVel(vec3<Scalar> U_t)
    {
    v_tan = U_t;
    }

    // // Accept optional velocity [RHEOINF]
    // /* \param  v_norm velocity
    //    from AnisoPotentialPair.h
    // */
    // HOSTDEVICE void setVel(vec3<Scalar> U_n)
    // {
    // v_norm = U_n;
    // }

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
        /*
        vec3<Scalar> rvec(dr);
        Scalar rsq = dot(rvec, rvec);
        Scalar r = fast::sqrt(rsq);
        //std::cout << "r..." << r << std::endl;
        Scalar radsum = 0.5 * (diameter_i + diameter_j); //~ Add radsum from passed diameters
        //std::cout << "radsum..." << radsum << std::endl;
        Scalar h_ij = r - radsum;
        //std::cout << "h_ij..." << h_ij << std::endl;

        //Scalar rcut = fast::sqrt(rcutsq);


        // if the particles do not interact
        if (rsq > rcutsq)
           {
           //std::cout << "Frix_false..." << std::endl;
           return false;
           }

        //Scalar morse_range = range;
        Scalar morse_range = 3/alpha;
        Scalar frix_cut = 0.5*morse_range;
        
        Scalar overlap = h_ij - frix_cut;

        vec3<Scalar> f;
        //f.x = 0;
        //f.y = 0;
        //f.z = 0;
        vec3<Scalar> t_i;
        //t_i.x = 0;
        //t_i.y = 0;
        //t_i.z = 0;
        vec3<Scalar> t_j;
        //t_j.x = 0;
        //t_j.y = 0;
        //t_j.z = 0;
        
        //std::cout << "overlap: " << overlap << std::endl;

        if (overlap < Scalar(0.0)) 
            {
            // Compute Friction
            // Scalar overlap = 1.0 - (h_ij / frix_cut); // 0 at h_ij=r_cut, 1 at colloid-colloid contact (prevents force blow up)
            

            vec3<Scalar> F_normal = - kn * overlap * rvec / r;  //repulsive since overlap = -ve
            // vec3<Scalar> F_normal = kn * overlap * rvec / r; //suggested: but + = attractive ; j->i
            
            //Scalar eta_n  = Scalar(0); 
            // vec3<Scalar> Fn_damp = eta_n * v_norm;//[AJ2024]

            vec3<Scalar> Ft_damp = -eta_t * v_tan;
            // vec3<Scalar> Ft_damp = eta_t * v_tan;
            
            Scalar F_normal_mag = fast::sqrt(dot(F_normal,F_normal));
            //Scalar Ft_normal_mag = mu_t * F_normal_mag;
            Scalar Ft_damp_mag =  fast::sqrt(dot(Ft_damp,Ft_damp));

            vec3<Scalar> t = normalize(v_tan);
            vec3<Scalar> friction_F = - std::min(mu_t*F_normal_mag,Ft_damp_mag) * t;
            //vec3<Scalar> friction_F = - Ft_damp_mag * t;
            // vec3<Scalar> friction_F = - mu_t*F_normal_mag * t;

            
            //std::cout << "Ft_damp..." << Ft_damp_mag << std::endl;
            //std::cout << "Ft_frix..." << Ft_normal_mag << std::endl;
            
            
            // std::cout << "Torque_x..." << t_i.x << std::endl;
            // std::cout << "Torque_x..." << t_i.x << std::endl;

            //f+= friction_F + Fn_damp + F_normal;
            //f+= friction_F + F_normal;
            f+=friction_F;

            //std::cout << "v_tan: " << v_tan.x << "," << v_tan.y << "," << v_tan.z <<  std::endl;
            //std::cout << "norm(v_tan): " << t.x << "," << t.y << "," << t.z <<  std::endl;
            //std::cout << "friction_F: " << friction_F.x << "," << friction_F.y << "," << friction_F.z <<  std::endl;
            t_i -= cross(0.5*(diameter_i+overlap)*rvec/r,friction_F);
            t_j -= cross(0.5*(diameter_j+overlap)*rvec/r,friction_F); 
                
            }

        //std::cout << "diameter_i: " << diameter_i << ", diameter_j: " << diameter_j << std::endl;
        //std::cout << "r: " << r << std::endl;
        //std::cout << "rvec: " << rvec.x << "," << rvec.y << "," << rvec.z << std::endl;
        //std::cout << "force: " << f.x << "," << f.y << "," << f.z <<  std::endl;
        //std::cout << "torque_i: " << t_i.x << "," << t_i.y << "," << t_i.z <<  std::endl;
        //std::cout << "torque_j: " << t_j.x << "," << t_j.y << "," << t_j.z <<  std::endl;
        //std::cout << "----" << std::endl;

            
        force = vec_to_scalar3(f);
        torque_i = vec_to_scalar3(t_i);
        torque_j = vec_to_scalar3(t_j);
        */

        vec3<Scalar> f_g;
        f_g.y = -F_G

        force = vec_to_scalar3(f_g)
        h = yposition 
        pair_eng += F_G*h //mgh 

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
        return "gravity";
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
    //Scalar range;
    //Scalar r0;     //!< Offset, i.e., position of the potential minimum
    Scalar kn;
    //Scalar eta_n;
    Scalar eta_t;
    Scalar mu_t;
    vec3<Scalar> v_tan;
    // vec3<Scalar> v_norm;        //[AJ2024]
    // const param_type &params;   //!< The pair potential parameters
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_MORSE_GRAVITY_H__
