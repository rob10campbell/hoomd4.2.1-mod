// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef __ANISO_POTENTIAL_PAIR_H__
#define __ANISO_POTENTIAL_PAIR_H__

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

#include "NeighborList.h"
#include "hoomd/ForceCompute.h"

#include "hoomd/ManagedArray.h"
#include "hoomd/VectorMath.h"

////~ angular rigidity [RHEOINF]
#include "Lifetime.h" //~ many-body neighbors [RHEOINF]

#include <pybind11/numpy.h>
//#include <pybind11/pybind11.h>
#include <stdexcept>

//#include "NeighborList.h"
//#include "hoomd/ForceCompute.h"
#include "hoomd/GlobalArray.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/Index1D.h"
#include "hoomd/managed_allocator.h"

//~ access for many-body neighbors [RHEOINF]
//#include "hoomd/HOOMDMPI.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
//~

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif
////~

//~ update some parameters every timestep (regardless of the neighborlist buffer) [RHEOINF]
#include "hoomd/HOOMDMPI.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#endif
//~

/*! \file AnisoPotentialPair.h
    \brief Defines the template class for anisotropic pair potentials
    \details The heart of the code that computes anisotropic pair potentials is in this file.
    \note This header cannot be compiled by nvcc
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace md
    {
//! Template class for computing pair potentials
/*! <b>Overview:</b>
    AnisoPotentialPair computes standard pair potentials (and forces) between all particle pairs in
   the simulation. It employs the use of a neighbor list to limit the number of computations done to
   only those particles with the cutoff radius of each other. The computation of the actual V(r) is
   not performed directly by this class, but by an aniso_evaluator class (e.g. EvaluatorPairLJ)
   which is passed in as a template parameter so the computations are performed as efficiently as
   possible.

    AnisoPotentialPair handles most of the gory internal details common to all standard pair
   potentials.
     - A cutoff radius to be specified per particle type pair
     - The energy can be globally shifted to 0 at the cutoff
     - Per type pair parameters are stored and a set method is provided
     - And all the details about looping through the particles, computing dr, computing the virial,
   etc. are handled

    \note XPLOR switching is not supported

    <b>Implementation details</b>

    rcutsq and the params are stored per particle type pair. It wastes a little bit of space, but
   benchmarks show that storing the symmetric type pairs and indexing with Index2D is faster than
   not storing redundant pairs and indexing with Index2DUpperTriangular. All of these values are
   stored in GlobalArray for easy access on the GPU by a derived class. The type of the parameters
   is defined by \a param_type in the potential aniso_evaluator class passed in. See the appropriate
   documentation for the aniso_evaluator for the definition of each element of the parameters.
*/

template<class aniso_evaluator> class AnisoPotentialPair : public ForceCompute
    {
    public:
    //! Param type from aniso_evaluator
    typedef typename aniso_evaluator::param_type param_type;

    //! Shape param type from aniso_evaluator
    typedef typename aniso_evaluator::shape_type shape_type;

    //! Construct the pair potential
    AnisoPotentialPair(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<NeighborList> nlist, //); //~ add angular rigidity params [RHEOINF]
                       Scalar K=0.0, Scalar w=0.0, Scalar theta_bar=0.0); //~ add angular rigidity params [RHEOINF]
    //! Destructor
    virtual ~AnisoPotentialPair();

    //! Set the pair parameters for a single type pair
    virtual void setParams(unsigned int typ1, unsigned int typ2, const param_type& param);

    virtual void setParamsPython(pybind11::tuple typ, pybind11::object params);

    /// Get params for a single type pair using a tuple of strings
    virtual pybind11::object getParamsPython(pybind11::tuple typ);

    //! Set the rcut for a single type pair
    virtual void setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut);

    /// Get the r_cut for a single type pair
    Scalar getRCut(pybind11::tuple types);

    /// Set the rcut for a single type pair using a tuple of strings
    virtual void setRCutPython(pybind11::tuple types, Scalar r_cut);

    /// Validate that types are within Ntypes
    virtual void validateTypes(unsigned int typ1, unsigned int typ2, std::string action);

    //! Set the shape parameters for a single type
    virtual void setShape(unsigned int typ, const shape_type& shape_param);

    virtual pybind11::object getShapePython(std::string typ);

    //! Set the shape parameters for a single type through Python
    virtual void setShapePython(std::string typ, const pybind11::object shape_param);

    std::vector<std::string> getTypeShapeMapping(
        const std::vector<param_type, hoomd::detail::managed_allocator<param_type>>& params,
        const std::vector<shape_type, hoomd::detail::managed_allocator<shape_type>>& shape_params)
        const
        {
        std::vector<std::string> type_shape_mapping(m_pdata->getNTypes());
        Scalar4 q = make_scalar4(1, 0, 0, 0);
        Scalar3 dr = make_scalar3(0, 0, 0);
        Scalar rcut = Scalar(0.0);
        for (unsigned int i = 0; i < type_shape_mapping.size(); i++)
            {
            aniso_evaluator evaluator(dr, q, q, rcut, m_params[m_typpair_idx(i, i)]);
            if (aniso_evaluator::needsShape())
                {
                evaluator.setShape(&shape_params[i], &shape_params[i]);
                }
            type_shape_mapping[i] = evaluator.getShapeSpec();
            }
        return type_shape_mapping;
        }

    pybind11::list getTypeShapesPy()
        {
        std::vector<std::string> type_shape_mapping
            = this->getTypeShapeMapping(m_params, m_shape_params);
        pybind11::list type_shapes;
        for (unsigned int i = 0; i < type_shape_mapping.size(); i++)
            type_shapes.append(type_shape_mapping[i]);
        return type_shapes;
        }

    //! Shifting modes that can be applied to the energy
    enum energyShiftMode
        {
        no_shift = 0,
        shift,
        };

    //! Set the mode to use for shifting the energy
    void setShiftMode(energyShiftMode mode)
        {
        m_shift_mode = mode;
        }

    void setShiftModePython(std::string mode)
        {
        if (mode == "none")
            {
            m_shift_mode = no_shift;
            }
        else if (mode == "shift")
            {
            m_shift_mode = shift;
            }
        else
            {
            throw std::runtime_error("Invalid energy shift mode.");
            }
        }

    /// Get the mod eused for the energy shifting
    std::string getShiftMode()
        {
        switch (m_shift_mode)
            {
        case no_shift:
            return "none";
        case shift:
            return "shift";
        default:
            return "";
            }
        }

    virtual void notifyDetach()
        {
        if (m_attached)
            {
            m_nlist->removeRCutMatrix(m_r_cut_nlist);
            }
        m_attached = false;
        }

    //~! Set angular rigidity params value [RHEOINF]
    void setK(Scalar K) {
        m_K = K;
    }
    Scalar getK() const {
        return m_K;
    }
    void setW(Scalar w) {
        m_w = w;
    }
    Scalar getW() const {
        return m_w;
    }
    void setTheta(Scalar theta_bar) {
        m_theta_bar = theta_bar;
    }
    Scalar getTheta() const {
        return m_theta_bar;
    }
    //~

    std::shared_ptr<Lifetime> LTIME; //~ access LTIME for many-body neighbors [RHEOINF]

#ifdef ENABLE_MPI
    //! Get ghost particle fields requested by this pair potential
    virtual CommFlags getRequestedCommFlags(uint64_t timestep);

    void updateGhostsIfNeeded(uint64_t timestep); //~ update some params every timestep regardless of buffer [RHEOINF]
#endif

    //! Returns true because we compute the torque
    virtual bool isAnisotropic()
        {
        return true;
        }

    /// Start autotuning kernel launch parameters
    virtual void startAutotuning()
        {
        ForceCompute::startAutotuning();

        // Start autotuning the neighbor list.
        m_nlist->startAutotuning();
        }

    /// Check if autotuning is complete.
    virtual bool isAutotuningComplete()
        {
        bool result = ForceCompute::isAutotuningComplete();
        result = result && m_nlist->isAutotuningComplete();
        return result;
        }

    protected:
    std::shared_ptr<NeighborList> m_nlist; //!< The neighborlist to use for the computation
    energyShiftMode m_shift_mode; //!< Store the mode with which to handle the energy shift at r_cut
    Scalar m_K;                   //!< angular rigidity K value [RHEOINF]
    Scalar m_w;                   //!< angular rigidity w value [RHEOINF]
    Scalar m_theta_bar;           //!< angular rigidity reference angle [RHEOINF]
    Index2D m_typpair_idx;        //!< Helper class for indexing per type pair arrays
    GlobalArray<Scalar> m_rcutsq; //!< Cutoff radius squared per type pair
    std::vector<param_type, hoomd::detail::managed_allocator<param_type>>
        m_params; //!< Pair parameters per type pair
    std::vector<shape_type, hoomd::detail::managed_allocator<shape_type>>
        m_shape_params; //!< Shape paramters per type

    /// Track whether we have attached to the Simulation object
    bool m_attached = true;

    /// r_cut (not squared) given to the neighbor list
    std::shared_ptr<GlobalArray<Scalar>> m_r_cut_nlist;

//~ update some params every timestep regardless of buffer [RHEOINF]
#ifdef ENABLE_MPI
   /// The system's communicator.
   std::shared_ptr<Communicator> m_comm;
#endif
//~

    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);
    };

/*! \param sysdef System to compute forces on
    \param nlist Neighborlist to use for computing the forces
*/
template<class aniso_evaluator>
AnisoPotentialPair<aniso_evaluator>::AnisoPotentialPair(std::shared_ptr<SystemDefinition> sysdef,
                                                        std::shared_ptr<NeighborList> nlist, //) //~ add angular rigidity params[RHEOINF]
                                                        Scalar K, Scalar w, Scalar theta_bar) //~ add angular rigidity params [RHEOINF]
    : ForceCompute(sysdef), m_nlist(nlist), m_shift_mode(no_shift),
      m_K(K), m_w(w), m_theta_bar(theta_bar), //~ add angular rigidity params [RHEOINF]
      m_typpair_idx(m_pdata->getNTypes())
    {
    m_exec_conf->msg->notice(5) << "Constructing AnisoPotentialPair<" << aniso_evaluator::getName()
                                << ">" << std::endl;
    assert(m_pdata);
    assert(m_nlist);

    //~ access Lifetime for indexing if needed [RHEOINF]
    if(m_K != 0.0)
        {
        LTIME = std::shared_ptr<Lifetime>(new Lifetime(sysdef));
        }
    //~

    GlobalArray<Scalar> rcutsq(m_typpair_idx.getNumElements(), m_exec_conf);
    m_rcutsq.swap(rcutsq);
    GlobalArray<Scalar> ronsq(m_typpair_idx.getNumElements(), m_exec_conf);
    std::vector<param_type, hoomd::detail::managed_allocator<param_type>> params(
        static_cast<size_t>(m_typpair_idx.getNumElements()),
        param_type(),
        hoomd::detail::managed_allocator<param_type>(m_exec_conf->isCUDAEnabled()));
    m_params.swap(params);

    std::vector<shape_type, hoomd::detail::managed_allocator<shape_type>> shape_params(
        static_cast<size_t>(m_pdata->getNTypes()),
        shape_type(),
        hoomd::detail::managed_allocator<shape_type>(m_exec_conf->isCUDAEnabled()));
    m_shape_params.swap(shape_params);

    m_r_cut_nlist
        = std::make_shared<GlobalArray<Scalar>>(m_typpair_idx.getNumElements(), m_exec_conf);
    nlist->addRCutMatrix(m_r_cut_nlist);

#if defined(ENABLE_HIP) && defined(__HIP_PLATFORM_NVCC__)
    if (m_exec_conf->isCUDAEnabled())
        {
        cudaMemAdvise(m_params.data(),
                      m_params.size() * sizeof(param_type),
                      cudaMemAdviseSetReadMostly,
                      0);
        cudaMemAdvise(m_shape_params.data(),
                      m_shape_params.size() * sizeof(shape_type),
                      cudaMemAdviseSetReadMostly,
                      0);

        // prefetch
        auto& gpu_map = m_exec_conf->getGPUIds();

        for (unsigned int idev = 0; idev < m_exec_conf->getNumActiveGPUs(); ++idev)
            {
            // prefetch data on all GPUs
            cudaMemPrefetchAsync(m_params.data(),
                                 sizeof(param_type) * m_params.size(),
                                 gpu_map[idev]);
            cudaMemPrefetchAsync(m_shape_params.data(),
                                 sizeof(shape_type) * m_shape_params.size(),
                                 gpu_map[idev]);
            }

        if (m_exec_conf->allConcurrentManagedAccess())
            {
            cudaMemAdvise(m_rcutsq.get(),
                          m_rcutsq.getNumElements() * sizeof(Scalar),
                          cudaMemAdviseSetReadMostly,
                          0);

            for (unsigned int idev = 0; idev < m_exec_conf->getNumActiveGPUs(); ++idev)
                {
                // prefetch data on all GPUs
                cudaMemPrefetchAsync(m_rcutsq.get(),
                                     sizeof(Scalar) * m_rcutsq.getNumElements(),
                                     gpu_map[idev]);
                }
            }
        }
#endif

//~ update some params every timestep regardless of neighborlist buffer [RHEOINF]
#ifdef ENABLE_MPI
    if (m_sysdef->isDomainDecomposed())
        {
        auto comm_weak = m_sysdef->getCommunicator();
        assert(comm_weak.lock());
        m_comm = comm_weak.lock();
        }
#endif
//~
    }

template<class aniso_evaluator> AnisoPotentialPair<aniso_evaluator>::~AnisoPotentialPair()
    {
    m_exec_conf->msg->notice(5) << "Destroying AnisoPotentialPair<" << aniso_evaluator::getName()
                                << ">" << std::endl;

    if (m_attached)
        {
        m_nlist->removeRCutMatrix(m_r_cut_nlist);
        }
    }

/*! \param typ1 First type index in the pair
    \param typ2 Second type index in the pair
    \param param Parameter to set
    \note When setting the value for (\a typ1, \a typ2), the parameter for (\a typ2, \a typ1) is
   automatically set.
*/
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setParams(unsigned int typ1,
                                                    unsigned int typ2,
                                                    const param_type& param)
    {
    validateTypes(typ1, typ2, "setting params");
    m_params[m_typpair_idx(typ1, typ2)] = param;
    m_params[m_typpair_idx(typ2, typ1)] = param;
    }

template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setParamsPython(pybind11::tuple typ,
                                                          pybind11::object params)
    {
    auto typ1 = m_pdata->getTypeByName(typ[0].cast<std::string>());
    auto typ2 = m_pdata->getTypeByName(typ[1].cast<std::string>());
    setParams(typ1, typ2, param_type(params, m_exec_conf->isCUDAEnabled()));
    }

template<class aniso_evaluator>
pybind11::object AnisoPotentialPair<aniso_evaluator>::getParamsPython(pybind11::tuple typ)
    {
    auto typ1 = m_pdata->getTypeByName(typ[0].cast<std::string>());
    auto typ2 = m_pdata->getTypeByName(typ[1].cast<std::string>());
    validateTypes(typ1, typ2, "getting params");

    return m_params[m_typpair_idx(typ1, typ2)].toPython();
    }

template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::validateTypes(unsigned int typ1,
                                                        unsigned int typ2,
                                                        std::string action)
    {
    // TODO change logic to just throw an exception
    auto n_types = this->m_pdata->getNTypes();
    if (typ1 >= n_types || typ2 >= n_types)
        {
        throw std::runtime_error("Error in" + action + " for pair potential. Invalid type");
        }
    }

/*! \param typ The type index.
    \param param Shape parameter to set
          set.
*/
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setShape(unsigned int typ, const shape_type& shape_param)
    {
    if (typ >= m_pdata->getNTypes())
        {
        throw std::runtime_error("Error setting shape parameters in AnisoPotentialPair");
        }

    m_shape_params[typ] = shape_param;
    }

/*! \param typ The type index.
    \param param Shape parameter to set
          set.
*/
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setShapePython(std::string typ,
                                                         pybind11::object shape_param)
    {
    auto typ_ = m_pdata->getTypeByName(typ);
    setShape(typ_, shape_type(shape_param, m_exec_conf->isCUDAEnabled()));
    }

/*! \param typ The type index.
    \param param Shape parameter to set
          set.
*/
template<class aniso_evaluator>
pybind11::object AnisoPotentialPair<aniso_evaluator>::getShapePython(std::string typ)
    {
    auto typ_ = m_pdata->getTypeByName(typ);
    if (typ_ >= m_pdata->getNTypes())
        {
        throw std::runtime_error("Error getting shape parameters in AnisoPotentialPair");
        }

    return m_shape_params[typ_].toPython();
    }

/*! \param typ1 First type index in the pair
    \param typ2 Second type index in the pair
    \param rcut Cutoff radius to set
    \note When setting the value for (\a typ1, \a typ2), the parameter for (\a typ2, \a typ1) is
   automatically set.
*/
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut)
    {
    validateTypes(typ1, typ2, "setting r_cut");
        {
        // store r_cut**2 for use internally
        ArrayHandle<Scalar> h_rcutsq(m_rcutsq, access_location::host, access_mode::readwrite);
        h_rcutsq.data[m_typpair_idx(typ1, typ2)] = rcut * rcut;
        h_rcutsq.data[m_typpair_idx(typ2, typ1)] = rcut * rcut;

        // store r_cut unmodified for so the neighbor list knows what particles to include
        ArrayHandle<Scalar> h_r_cut_nlist(*m_r_cut_nlist,
                                          access_location::host,
                                          access_mode::readwrite);
        h_r_cut_nlist.data[m_typpair_idx(typ1, typ2)] = rcut;
        h_r_cut_nlist.data[m_typpair_idx(typ2, typ1)] = rcut;
        }

    // notify the neighbor list that we have changed r_cut values
    m_nlist->notifyRCutMatrixChange();
    }

template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::setRCutPython(pybind11::tuple types, Scalar r_cut)
    {
    auto typ1 = m_pdata->getTypeByName(types[0].cast<std::string>());
    auto typ2 = m_pdata->getTypeByName(types[1].cast<std::string>());
    setRcut(typ1, typ2, r_cut);
    }

template<class aniso_evaluator>
Scalar AnisoPotentialPair<aniso_evaluator>::getRCut(pybind11::tuple types)
    {
    auto typ1 = m_pdata->getTypeByName(types[0].cast<std::string>());
    auto typ2 = m_pdata->getTypeByName(types[1].cast<std::string>());
    validateTypes(typ1, typ2, "getting r_cut.");
    ArrayHandle<Scalar> h_rcutsq(m_rcutsq, access_location::host, access_mode::read);
    return sqrt(h_rcutsq.data[m_typpair_idx(typ1, typ2)]);
    }

/*! \post The pair forces are computed for the given timestep. The neighborlist's compute method is
   called to ensure that it is up to date before proceeding.

    \param timestep specifies the current time step of the simulation
*/
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::computeForces(uint64_t timestep)
    {
    // start by updating the neighborlist
    m_nlist->compute(timestep);

    //~ update some params every timestep regardless of neighborlist buffer [RHEOINF]
    m_nlist->setStorageMode(NeighborList::full);
    //~

    // depending on the neighborlist settings, we can take advantage of newton's third law
    // to reduce computations at the cost of memory access complexity: set that flag now
    bool third_law = m_nlist->getStorageMode() == NeighborList::half;

    // access the neighbor list, particle data, and system box
    ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(),
                                        access_location::host,
                                        access_mode::read);
    ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(),
                                      access_location::host,
                                      access_mode::read);
    ArrayHandle<size_t> h_head_list(m_nlist->getHeadList(),
                                    access_location::host,
                                    access_mode::read);

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    //~ add diameter [RHEOINF]
    ArrayHandle<Scalar> h_diameter(m_pdata->getDiameters(),
                                   access_location::host,
                                   access_mode::read);
    //~
    ArrayHandle<Scalar> h_charge(m_pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_orientation(m_pdata->getOrientationArray(),
                                       access_location::host,
                                       access_mode::read);
    //~ add velocity [RHEOINF]
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(),
                               access_location::host,
                               access_mode::read);
    //~
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);

    //~ add angular momentul and intertia [RHEOINF]
    ArrayHandle<Scalar4> h_angmom(m_pdata->getAngularMomentumArray(),
                                  access_location::host,
                                  access_mode::read);
    ArrayHandle<Scalar3> h_inertia(m_pdata->getMomentsOfInertiaArray(),
                                   access_location::host,
                                   access_mode::read);
    //~

    // force arrays
    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar4> h_torque(m_torque, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar> h_virial(m_virial, access_location::host, access_mode::overwrite);

    //~ update some params every timestep (regardless of neighborlist buffer) [RHEOINF]
    #ifdef ENABLE_MPI
        updateGhostsIfNeeded(timestep);
    #endif
    //~

    const BoxDim box = m_pdata->getBox();
    ArrayHandle<Scalar> h_rcutsq(m_rcutsq, access_location::host, access_mode::read);
        {
        // need to start from a zero force, energy and virial
        memset(&h_force.data[0], 0, sizeof(Scalar4) * m_pdata->getN());
        memset(&h_torque.data[0], 0, sizeof(Scalar4) * m_pdata->getN());
        memset(&h_virial.data[0], 0, sizeof(Scalar) * m_virial.getNumElements());

        PDataFlags flags = this->m_pdata->getFlags();
        bool compute_virial = flags[pdata_flag::pressure_tensor];

        //~ get multi-body neighbors [RHEOINF]
        //S
        // NOTE: h_diameter is already included above
        ArrayHandle<unsigned int> h_tag(m_pdata->getTags(),
                                        access_location::host,
                                        access_mode::read);

        ArrayHandle<Scalar> h_current_neighbor_list(m_pdata->getParticleNList(),
                                                        access_location::host,
                                                        access_mode::readwrite);
        assert(h_current_neighbor_list.data);
        size_t p_neighbor_pitch = m_pdata->getParticleNList().getPitch();
        //std::cout<<"new 2"<<std::endl;

        //~ update ghost neighbors [RHEOINF]
        int tot_particles = (int)(m_pdata->getN() + m_pdata->getNGhosts());
        std::vector<Scalar> h_previous_neighbor_list(20 * tot_particles);
        #ifdef ENABLE_MPI
            if (m_K != 0.0 ){
            updateGhostsIfNeeded(timestep);}
        #endif

        if (m_K != 0.0 ){
        
            for ( int i = 0; i < tot_particles; i++) {
                for (size_t j = 0; j < 20; j++) { // 20 is the count of neighbors used
                    h_previous_neighbor_list[j * (tot_particles) + i] = h_current_neighbor_list.data[j * p_neighbor_pitch + i];
                }
            }
    
        }

        // start the connected neighbors with -2 
        for (int i = 0; i <  tot_particles ; ++i) {
            for(int j = 0; j <  20 ; ++j){ 
            h_current_neighbor_list.data[j * p_neighbor_pitch + i] = -2 ;
            }
        }
        // start the neighbors with tag of the particle so that the particle can be found
        for (int i = 0; i < (int)m_pdata->getN(); i++){
            h_current_neighbor_list.data[i] = h_tag.data[i];
    
        }
        //F
        //~

        size_t idx_pi = -1; //~ set the neighbor index for pi [RHEOINF]

        // for each particle
        for (int i = 0; i < (int)m_pdata->getN(); i++)
            {
            // access the particle's position and type (MEM TRANSFER: 4 scalars)
            Scalar3 pi = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
            unsigned int typei = __scalar_as_int(h_pos.data[i].w);
            Scalar4 quat_i = h_orientation.data[i];

            //~ update many-body neighbors [RHEOINF]
            if (m_K != 0.0){
                //find the particle in the previous time step particle neighbor list
                if(typei == 0 ){
                    //bool wasnt_found = true;
                    for(int p_i = 0; p_i < tot_particles; p_i++){
                        if (h_previous_neighbor_list[p_i] == h_tag.data[i]) {
                            idx_pi = p_i;
                            //wasnt_found = false;
                            break;
                        }else{idx_pi = -3;}
                    }
                }
            }
            //~

            //~ angular momentum / quaternion updates [RHEOINF]
            Scalar3 vi = make_scalar3(h_vel.data[i].x, h_vel.data[i].y, h_vel.data[i].z);
            quat<Scalar> p(h_angmom.data[i]);
            quat<Scalar> q(h_orientation.data[i]);
            vec3<Scalar> I(h_inertia.data[i]);
            vec3<Scalar> omega_i = (Scalar(1. / 2.) * conj(q) * p).v / I;
            //std::cout << typei << "," << vi.x << "," << vi.y << "," << vi.z << std::endl;
            //~

            // sanity check
            assert(typei < m_pdata->getNTypes());

            // access charge (if needed)
            //~ add diameter [RHEOINF]
            Scalar di = Scalar(0.0);
            //~
            Scalar qi = Scalar(0.0);
            //~ add diameter [RHEOINF]
            if (aniso_evaluator::needsDiameter())
                di = h_diameter.data[i];
            //~
            if (aniso_evaluator::needsCharge())
                qi = h_charge.data[i];

            // initialize current particle force, torque, potential energy, and virial to 0
            Scalar fxi = Scalar(0.0);
            Scalar fyi = Scalar(0.0);
            Scalar fzi = Scalar(0.0);
            Scalar txi = Scalar(0.0);
            Scalar tyi = Scalar(0.0);
            Scalar tzi = Scalar(0.0);
            Scalar pei = Scalar(0.0);
            Scalar virialxxi = 0.0;
            Scalar virialxyi = 0.0;
            Scalar virialxzi = 0.0;
            Scalar virialyyi = 0.0;
            Scalar virialyzi = 0.0;
            Scalar virialzzi = 0.0;

            // loop over all of the neighbors of this particle
            const size_t myHead = h_head_list.data[i];
            const unsigned int size = (unsigned int)h_n_neigh.data[i];
            for (unsigned int k = 0; k < size; k++)
                {
                // access the index of this neighbor (MEM TRANSFER: 1 scalar)
                unsigned int j = h_nlist.data[myHead + k];
                assert(j < m_pdata->getN() + m_pdata->getNGhosts());

                // calculate dr_ji (MEM TRANSFER: 3 scalars / FLOPS: 3)
                Scalar3 pj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
                Scalar3 dx = pi - pj;
                Scalar4 quat_j = h_orientation.data[j];

                // access the type of the neighbor particle (MEM TRANSFER: 1 scalar)
                unsigned int typej = __scalar_as_int(h_pos.data[j].w);
                assert(typej < m_pdata->getNTypes());

                // access charge (if needed)
                //~ add diameter [RHEOINF]
                Scalar dj = Scalar(0.0);
                //~
                Scalar qj = Scalar(0.0);
                //~ add diameter [RHEOINF]
                if (aniso_evaluator::needsDiameter())
                    dj = h_diameter.data[j];
                //~
                if (aniso_evaluator::needsCharge())
                    qj = h_charge.data[j];

                // apply periodic boundary conditions
                dx = box.minImage(dx);

                //~ calculate r_ij squared (FLOPS: 5) [RHEOINF]
                Scalar rsq = dot(dx, dx);
                //~

                //~ calculate the center-center distance equal to particle-particle contact (AKA r0) [RHEOINF]
                //Scalar radcontact = Scalar(0.5) * (h_diameter.data[i] + h_diameter.data[j]);
                //~

                //~ angular momentum / quaternion updates [RHEOINF]
                vec3<Scalar> dx_vec(dx);
                vec3<Scalar> n_normal = normalize(dx_vec);
                Scalar3 vj = make_scalar3(h_vel.data[j].x, h_vel.data[j].y, h_vel.data[j].z);
                vec3<Scalar> vij(vi-vj);
                quat<Scalar> p(h_angmom.data[j]);
                quat<Scalar> q(h_orientation.data[j]);
                vec3<Scalar> I(h_inertia.data[j]);
                vec3<Scalar> omega_j = (Scalar(1. / 2.) * conj(q) * p).v / I;
                vec3<Scalar> Vt = vij - dot(vij,n_normal) * n_normal - Scalar(0.5) * cross((di*omega_i+dj*omega_j),n_normal);
                //~

                // get parameters for this type pair
                unsigned int typpair_idx = m_typpair_idx(typei, typej);
                const param_type& param = m_params[typpair_idx];
                Scalar rcutsq = h_rcutsq.data[typpair_idx];

                // design specifies that energies are shifted if
                // shift mode is set to shift
                bool energy_shift = false;
                if (m_shift_mode == shift)
                    energy_shift = true;

                // compute the force and potential energy
                Scalar3 force = make_scalar3(0.0, 0.0, 0.0);
                Scalar3 torque_i = make_scalar3(0.0, 0.0, 0.0);
                Scalar3 torque_j = make_scalar3(0.0, 0.0, 0.0);

                Scalar pair_eng = Scalar(0.0);

                aniso_evaluator eval(dx, quat_i, quat_j, rcutsq, param);

                //~ add diameter [RHEOINF]
                if (aniso_evaluator::needsDiameter())
                    eval.setDiameter(di, dj);
                //~
                if (aniso_evaluator::needsCharge())
                    eval.setCharge(qi, qj);
                if (aniso_evaluator::needsShape())
                    eval.setShape(&m_shape_params[typei], &m_shape_params[typej]);
                if (aniso_evaluator::needsTags())
                    eval.setTags(h_tag.data[i], h_tag.data[j]);
                //~ add pair typeids for polydispersity [RHEOINF]
                if (aniso_evaluator::needsTypes())
                    eval.setTypes(typei, typej);
                //~
                //~ add velocity [RHEOINF]
                if (aniso_evaluator::needsVel())
                    eval.setVel(Vt);
                //~


                //~ update multi-body angles [RHEOINF]
                //S
                unsigned int tagi = h_tag.data[i];
                unsigned int tagj = h_tag.data[j];

                if (m_K != 0.0){
                    if (typei == typej)
                    {
                        Scalar rsq_root = std::sqrt(rsq) - Scalar(2.0);
                        if (rsq_root < Scalar(2.0))
                        {
                            // Check and update neighbors for particle i
                            size_t idx_i = 0;
                            bool not_saved = true;
                            while (idx_i < 20){

                                // Check if tag[j] already exists
                                if (h_current_neighbor_list.data[idx_i * p_neighbor_pitch + i] == h_tag.data[j])
                                {
                                    not_saved = false;
                                    break;
                                }

                                // If current entry is empty, save tag[j]
                                if (h_current_neighbor_list.data[idx_i * p_neighbor_pitch + i] == -2)
                                {
                                    h_current_neighbor_list.data[idx_i * p_neighbor_pitch + i] = h_tag.data[j];
                                    break;
                                }

                                idx_i++;
                            }

                            if (not_saved)
                            {
                                // Check and update neighbors for particle j
                                size_t idx_j = 0;
                                while (idx_j < 20)
                                {
                                    // Check if tag[i] already exists
                                    if (h_current_neighbor_list.data[idx_j * p_neighbor_pitch + j] == h_tag.data[i])
                                    {
                                        break;
                                    }

                                    // If current entry is empty, save tag[i]
                                    if (h_current_neighbor_list.data[idx_j * p_neighbor_pitch + j] == -2)
                                    {
                                        h_current_neighbor_list.data[idx_j * p_neighbor_pitch + j] = h_tag.data[i];
                                        break;
                                    }

                                    idx_j++;
                                }
                            }


                                // save angles 
                                bool new_connection = true;
                                size_t is_idx = 0;
                                while ( is_idx < 20)
                                {
                                    if (h_previous_neighbor_list[is_idx * tot_particles + idx_pi] == h_tag.data[j])
                                    {
                                        new_connection = false;
                                        break;
                                    }
                                    is_idx++;
                                }

                                
                        }           
                     }                   
                 }               
                 //F
                 //~


                bool evaluated = eval.evaluate(force, pair_eng, energy_shift, torque_i, torque_j);

                if (evaluated)
                    {
                    Scalar3 force2 = Scalar(0.5) * force;

                    // add the force, potential energy and virial to the particle i
                    // (FLOPS: 8)
                    fxi += force.x;
                    fyi += force.y;
                    fzi += force.z;
                    txi += torque_i.x;
                    tyi += torque_i.y;
                    tzi += torque_i.z;
                    pei += pair_eng * Scalar(0.5);

                    if (compute_virial)
                        {
                        virialxxi += dx.x * force2.x;
                        virialxyi += dx.y * force2.x;
                        virialxzi += dx.z * force2.x;
                        virialyyi += dx.y * force2.y;
                        virialyzi += dx.z * force2.y;
                        virialzzi += dx.z * force2.z;
                        }

                    // add the force to particle j if we are using the third law (MEM TRANSFER: 10
                    // scalars / FLOPS: 8)
                    if (third_law)
                        {
                        h_force.data[j].x -= force.x;
                        h_force.data[j].y -= force.y;
                        h_force.data[j].z -= force.z;
                        h_torque.data[j].x += torque_j.x;
                        h_torque.data[j].y += torque_j.y;
                        h_torque.data[j].z += torque_j.z;
                        h_force.data[j].w += pair_eng * Scalar(0.5);
                        if (compute_virial)
                            {
                            h_virial.data[0 * m_virial_pitch + j] += dx.x * force2.x;
                            h_virial.data[1 * m_virial_pitch + j] += dx.y * force2.x;
                            h_virial.data[2 * m_virial_pitch + j] += dx.z * force2.x;
                            h_virial.data[3 * m_virial_pitch + j] += dx.y * force2.y;
                            h_virial.data[4 * m_virial_pitch + j] += dx.z * force2.y;
                            h_virial.data[5 * m_virial_pitch + j] += dx.z * force2.z;
                            }
                        }
                    }
                }

            // finally, increment the force, potential energy and virial for particle i
            h_force.data[i].x += fxi;
            h_force.data[i].y += fyi;
            h_force.data[i].z += fzi;
            h_torque.data[i].x += txi;
            h_torque.data[i].y += tyi;
            h_torque.data[i].z += tzi;
            h_force.data[i].w += pei;
            if (compute_virial)
                {
                h_virial.data[0 * m_virial_pitch + i] += virialxxi;
                h_virial.data[1 * m_virial_pitch + i] += virialxyi;
                h_virial.data[2 * m_virial_pitch + i] += virialxzi;
                h_virial.data[3 * m_virial_pitch + i] += virialyyi;
                h_virial.data[4 * m_virial_pitch + i] += virialyzi;
                h_virial.data[5 * m_virial_pitch + i] += virialzzi;
                }
            }

           //~ compute angular rigidity forces [RHEOINF]
           //S
if (m_K != 0.0){
    #ifdef ENABLE_MPI

        LTIME->updatebondtime(timestep); //~ get timestep from Liftime file [RHEOINF]

    #endif
    for (int i = 0; i < (int)m_pdata->getN(); i++){
        unsigned int typei = __scalar_as_int(h_pos.data[i].w);
        Scalar3 pi = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
        if (typei==0) {
                //Scalar kk =0;
                size_t nonzero_count = 0;
                unsigned int tagi = h_tag.data[i];
                for (size_t idx_i = 1; idx_i < 20 ; ++idx_i)
                {
                    // check if i has more than 1 neighbor
                    if (h_current_neighbor_list.data[idx_i* p_neighbor_pitch +i] != -2
                    && h_current_neighbor_list.data[idx_i* p_neighbor_pitch +i] != tagi)
                    {
                        nonzero_count++; //count the neighbors
                    }
                }
                if (nonzero_count > 1)
                {
                    #define SMALL Scalar(0.001)
                    // go over each angle between two different neighbors of i and apply angle force 
                    // find the first neighbor a
                    for (size_t idx_a = 1; idx_a < (nonzero_count); ++idx_a)
                    {
                        unsigned int taga = static_cast<unsigned int>(h_current_neighbor_list.data[idx_a* p_neighbor_pitch +i]);
                        Scalar3 pa;
                        int a = -1;
                        //bool f5=false;
                        for (int M = 0; M < (int)(m_pdata->getN()+ m_pdata->getNGhosts()); M++)
                        {
                            if (h_tag.data[M] == taga) {
                            // Access the position of the particle with the matching tag
                                pa = make_scalar3(h_pos.data[M].x, h_pos.data[M].y, h_pos.data[M].z);
                                a = M;
                                //f5=true;
                                break;
                                }
                        }

                        unsigned int typea = __scalar_as_int(h_pos.data[a].w);
                        // find the second neighbor b
                        for (size_t idx_b = (idx_a+1) ; idx_b < (nonzero_count+1); ++idx_b)
                        {
                            unsigned int tagb = static_cast<unsigned int>(h_current_neighbor_list.data[idx_b* p_neighbor_pitch +i]);
                            Scalar3 pb;
                            int b = -1;
                            //bool f6=false;
                            for (int N = 0; N < (int)(m_pdata->getN()+ m_pdata->getNGhosts()); N++)
                            {
                                if (h_tag.data[N] == tagb) {
                                    // Access the position of the particle with the matching tag
                                    pb = make_scalar3(h_pos.data[N].x, h_pos.data[N].y, h_pos.data[N].z);
                                    b = N;
                                    //f6 = true;
                                    break; 
                                    }

                            }
                            unsigned int typeb = __scalar_as_int(h_pos.data[b].w);



                            if (tagb != taga  && taga != tagi && tagb != tagi && typea==0 && typeb==0)
                            {
                                // Calculate angle index

                                unsigned int vari = tagi ;
                                unsigned int vara = taga ;
                                unsigned int varb = tagb ;

                                unsigned int var1 = vari;
                                unsigned int var3 = std::max({vara, varb});
                                unsigned int var2 = std::min({vara, varb});

                                unsigned int n = LTIME->num_solvent;
                                unsigned int current_angle_index = (var1 * n*(n-1)/2) + (2*var2*n - var2*var2 + 2*var3 - 3*var2 -2)/2;

                                Scalar3 ria = make_scalar3(pa.x - pi.x, pa.y - pi.y, pa.z - pi.z);
                                Scalar3 rib = make_scalar3(pb.x - pi.x, pb.y - pi.y, pb.z - pi.z);
                                ria = box.minImage(ria);
                                rib = box.minImage(rib);

                                // Calculate magnitudes of vectors
                                Scalar ria_mag = std::sqrt(ria.x * ria.x + ria.y * ria.y + ria.z * ria.z);
                                Scalar rib_mag = std::sqrt(rib.x * rib.x + rib.y * rib.y + rib.z * rib.z);

                                // Calculate dot product
                                Scalar dot_product = ria.x * rib.x + ria.y * rib.y + ria.z * rib.z;

                                // Calculate cosine of the angle 
                                Scalar current_cos_theta = dot_product / (ria_mag * rib_mag);


                                if (current_cos_theta > 1.0)
                                    current_cos_theta = 1.0;
                                if (current_cos_theta < -1.0)
                                    current_cos_theta = -1.0;

                                Scalar current_sin_theta = std::sqrt(1.0 - current_cos_theta * current_cos_theta);
                                if (current_sin_theta < SMALL)
                                    current_sin_theta = SMALL;

                                current_sin_theta = 1.0 / current_sin_theta;

                                // Calculate force magnitude for Harmonic equation
                                //Scalar m_theta_bar = 120;
                                Scalar eq_theta = m_theta_bar * M_PI / 180; //65.0 * M_PI / 180.0;
                                Scalar dth = acos(current_cos_theta) - eq_theta;
                                //Scalar m_K = 100;
                                Scalar tk = m_K * dth;

                                //Calculate force magnitude for K(thera - theta_0)^2
                                Scalar vab = -1.0 * tk * current_sin_theta;

                                Scalar a11 = vab * current_cos_theta / (ria_mag * ria_mag);
                                Scalar a12 = -vab / (ria_mag * rib_mag);
                                Scalar a22 = vab * current_cos_theta / (rib_mag * rib_mag);

                                // Calculate forces
                                Scalar3 fia = make_scalar3(a11 * ria.x + a12 * rib.x,
                                                            a11 * ria.y + a12 * rib.y,
                                                            a11 * ria.z + a12 * rib.z);

                                Scalar3 fib = make_scalar3(a22 * ria.x + a12 * rib.x,
                                                            a22 * ria.y + a12 * rib.y,
                                                            a22 * ria.z + a12 * rib.z);

                                // compute the energy, for each atom in the angle for K(thera - theta_0)^2
                                Scalar angle_eng = (tk * dth) * Scalar(1.0 / 6.0);

                                Scalar angle_virial[6];
                                angle_virial[0] = Scalar(1. / 3.) * (ria.x * fia.x + rib.x * fib.x);
                                angle_virial[1] = Scalar(1. / 3.) * (ria.y * fia.x + rib.y * fib.x);
                                angle_virial[2] = Scalar(1. / 3.) * (ria.z * fia.x + rib.z * fib.x);
                                angle_virial[3] = Scalar(1. / 3.) * (ria.y * fia.y + rib.y * fib.y);
                                angle_virial[4] = Scalar(1. / 3.) * (ria.z * fia.y + rib.z * fib.y);
                                angle_virial[5] = Scalar(1. / 3.) * (ria.z * fia.z + rib.z * fib.z);


                                // Update forces and virials for particle i, a and b
                                if (a < (int)(m_pdata->getN()+ m_pdata->getNGhosts())) {
                                    h_force.data[a].x += fia.x;
                                    h_force.data[a].y += fia.y;
                                    h_force.data[a].z += fia.z;
                                    h_force.data[a].w += angle_eng;
                                    for (int l = 0; l < 6; l++)
                                        h_virial.data[l * m_virial_pitch + a] += angle_virial[l];

                                }

                                if (i < (int)m_pdata->getN()) {
                                    h_force.data[i].x -= fia.x + fib.x;
                                    h_force.data[i].y -= fia.y + fib.y;
                                    h_force.data[i].z -= fia.z + fib.z;
                                    h_force.data[i].w += angle_eng;
                                    for (int l = 0; l < 6; l++)
                                        h_virial.data[l * m_virial_pitch + i] += angle_virial[l];

                                }

                                if (b < (int)(m_pdata->getN()+ m_pdata->getNGhosts())) {
                                    h_force.data[b].x += fib.x;
                                    h_force.data[b].y += fib.y;
                                    h_force.data[b].z += fib.z;
                                    h_force.data[b].w += angle_eng;
                                    for (int l = 0; l < 6; l++)
                                        h_virial.data[l * m_virial_pitch + b] += angle_virial[l];
                                }

                            }
                        }
                    }
                }
            }
        }
}
//F


        }
    }


//~ update s me parameters every timestep (regardless of neighborlist buffer) [RHEOINF]
#ifdef ENABLE_MPI
template<class aniso_evaluator>
void AnisoPotentialPair<aniso_evaluator>::updateGhostsIfNeeded(uint64_t timestep)
    {
    // Temporarily modify communication flags to ensure all required data is exchanged
    CommFlags old_flags = m_comm->getFlags();
    CommFlags new_flags = old_flags;
    new_flags[comm_flag::tag] = 1;
    new_flags[comm_flag::velocity] = 1;
    new_flags[comm_flag::orientation] = 1;
    if (aniso_evaluator::needsDiameter())
        new_flags[comm_flag::diameter] = 1;
    new_flags[comm_flag::net_torque] = 1;

    m_comm->setFlags(new_flags);

    // Force communication
    m_comm->migrateParticles();
    m_comm->exchangeGhosts();

    // Restore original flags
    m_comm->setFlags(old_flags);
    }
#endif
//~

#ifdef ENABLE_MPI
/*! \param timestep Current time step
 */
template<class aniso_evaluator>
CommFlags AnisoPotentialPair<aniso_evaluator>::getRequestedCommFlags(uint64_t timestep)
    {
    CommFlags flags = CommFlags(0);

    //~ we need velocity [RHEOINF]
    flags[comm_flag::velocity] = 1;
    //~

    // we need orientations for anisotropic ptls
    flags[comm_flag::orientation] = 1;

    if (aniso_evaluator::needsCharge())
        flags[comm_flag::charge] = 1;

    //~ add diameter [RHEOINF]
    if (aniso_evaluator::needsDiameter())
        flags[comm_flag::diameter] = 1;
    //~

    // with rigid bodies, include net torque
    flags[comm_flag::net_torque] = 1;

    flags |= ForceCompute::getRequestedCommFlags(timestep);

    return flags;
    }
#endif

namespace detail
    {
//! Export this pair potential to python
/*! \param name Name of the class in the exported python module
    \tparam T Evaluator type to export.
*/
template<class T> void export_AnisoPotentialPair(pybind11::module& m, const std::string& name)
    {
    pybind11::class_<AnisoPotentialPair<T>, ForceCompute, std::shared_ptr<AnisoPotentialPair<T>>>
        anisopotentialpair(m, name.c_str());
    anisopotentialpair
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, Scalar, Scalar, Scalar>()) //~ add Scalar for angular repulsion params [RHEOINF]
        .def("setParams", &AnisoPotentialPair<T>::setParamsPython)
        .def("getParams", &AnisoPotentialPair<T>::getParamsPython)
        .def("setShape", &AnisoPotentialPair<T>::setShapePython)
        .def("getShape", &AnisoPotentialPair<T>::getShapePython)
        .def("setRCut", &AnisoPotentialPair<T>::setRCutPython)
        .def("getRCut", &AnisoPotentialPair<T>::getRCut)
        .def_property("mode",
                      &AnisoPotentialPair<T>::getShiftMode,
                      &AnisoPotentialPair<T>::setShiftModePython)
        .def_property("K", &AnisoPotentialPair<T>::getK, &AnisoPotentialPair<T>::setK) //~ add K [RHEOINF]
        .def_property("w", &AnisoPotentialPair<T>::getW, &AnisoPotentialPair<T>::setW) //~ add w [RHEOINF]
        .def_property("theta_bar", &AnisoPotentialPair<T>::getTheta, &AnisoPotentialPair<T>::setTheta) //~ add theta_bar [RHEOINF]
        .def("getTypeShapesPy", &AnisoPotentialPair<T>::getTypeShapesPy);
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd

#endif // __ANISO_POTENTIAL_PAIR_H__
