// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "ConstraintAngleForceCompute.h"

#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>

using namespace std;

// SMALL a relatively small number
#define SMALL Scalar(0.001)

/*! \file ConstraintAngleForceCompute.cc
    \brief Contains code for the ConstraintAngleForceCompute class
*/

namespace hoomd
    {
namespace md
    {
/*! \param sysdef System to compute forces on
    \post Memory is allocated, and forces are zeroed.
*/
ConstraintAngleForceCompute::ConstraintAngleForceCompute(std::shared_ptr<SystemDefinition> sysdef,
							std::shared_ptr<NeighborList> nlist)
    : ForceCompute(sysdef), m_nlist(nlist), m_K(NULL), m_t_0(NULL)
    {
    m_exec_conf->msg->notice(5) << "Constructing ConstraintAngleForceCompute" << endl;

    assert(m_nlist);
    assert(m_pdata);

    //m_nlist->setStorageMode(NeighborList::full);
    // access the angle data for later use
    m_angle_data = m_sysdef->getAngleData();

    // check for some silly errors a user could make
    if (m_angle_data->getNTypes() == 0)
        {
        throw runtime_error("No angle types in the system.");
        }

    // allocate the parameters
    m_K = new Scalar[m_angle_data->getNTypes()];
    m_t_0 = new Scalar[m_angle_data->getNTypes()];
    }

ConstraintAngleForceCompute::~ConstraintAngleForceCompute()
    {
    m_exec_conf->msg->notice(5) << "Destroying ConstraintAngleForceCompute" << endl;

    delete[] m_K;
    delete[] m_t_0;
    m_K = NULL;
    m_t_0 = NULL;
    }

/*! \param type Type of the angle to set parameters for
    \param K Stiffness parameter for the force computation
    \param t_0 Equilibrium angle in radians for the force computation

    Sets parameters for the potential of a particular angle type
*/
void ConstraintAngleForceCompute::setParams(unsigned int type, Scalar K, Scalar t_0)
    {
    // make sure the type is valid
    if (type >= m_angle_data->getNTypes())
        {
        throw runtime_error("Invalid angle type.");
        }

    m_K[type] = K;
    m_t_0[type] = t_0;

    // check for some silly errors a user could make
    if (K <= 0)
        m_exec_conf->msg->warning() << "angle.constraint: specified K <= 0" << endl;
    if (t_0 <= 0)
        m_exec_conf->msg->warning() << "angle.constraint: specified t_0 <= 0" << endl;
    }

void ConstraintAngleForceCompute::setParamsPython(std::string type, pybind11::dict params)
    {
    auto typ = m_angle_data->getTypeByName(type);
    auto _params = angle_harmonic_params(params);
    setParams(typ, _params.k, _params.t_0);
    }

pybind11::dict ConstraintAngleForceCompute::getParams(std::string type)
    {
    auto typ = m_angle_data->getTypeByName(type);
    if (typ >= m_angle_data->getNTypes())
        {
        throw runtime_error("Invalid angle type.");
        }
    pybind11::dict params;
    params["k"] = m_K[typ];
    params["t0"] = m_t_0[typ];
    return params;
    }

/*! Actually perform the force computation
    \param timestep Current time step
 */
void ConstraintAngleForceCompute::computeForces(uint64_t timestep)
    {
    assert(m_pdata);
    m_nlist->compute(timestep);
    // access the particle data arrays
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

    ArrayHandle<Scalar> h_diameter(m_pdata->getDiameters(),
                                   access_location::host,
                                   access_mode::read);

    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar> h_virial(m_virial, access_location::host, access_mode::overwrite);
    size_t virial_pitch = m_virial.getPitch();

    // there are enough other checks on the input data: but it doesn't hurt to be safe
    assert(h_force.data);
    assert(h_virial.data);
    assert(h_pos.data);
    assert(h_rtag.data);

    // Zero data for force calculation.
    memset((void*)h_force.data, 0, sizeof(Scalar4) * m_force.getNumElements());
    memset((void*)h_virial.data, 0, sizeof(Scalar) * m_virial.getNumElements());

    // get a local copy of the simulation box too
    const BoxDim& box = m_pdata->getGlobalBox();

    ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(),
                                        access_location::host,
                                        access_mode::read);
    ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(),
                                      access_location::host,
                                      access_mode::read);
    ArrayHandle<size_t> h_head_list(m_nlist->getHeadList(),
                                    access_location::host,
                                    access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(),
                                    access_location::host,
                                    access_mode::read);

    std::vector<std::vector<unsigned int> > TempList((int)m_pdata->getN());
    unsigned int i_temp = 0;
    for (int i = 0; i < (int)m_pdata->getN(); i++)
        {
        unsigned int tagi = h_tag.data[i];
        Scalar3 pi = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
        const size_t myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        i_temp += size;
        for (unsigned int k = 0; k < size; k++)
            {
            unsigned int j = h_nlist.data[myHead + k];
            assert(j < m_pdata->getN() + m_pdata->getNGhosts());
            Scalar3 pj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
            Scalar3 dx = pi - pj;
            dx = box.minImage(dx);
            Scalar hij = fast::sqrt(dot(dx, dx))-0.5*(h_diameter.data[i]+h_diameter.data[j]);
            if(hij<0.08)
                {
                TempList[i].push_back(h_tag.data[j]);
                TempList[j].push_back(tagi);
                }
            }
        }

    for (int i = 0; i < (int)m_pdata->getN(); i++)
        {
        if(TempList[i].size()>=2)
            {
            sort(TempList[i].begin(), TempList[i].end());
            for (unsigned int j = 0; j < TempList[i].size()-1; j++)
                {
                for (unsigned int k = j+1; k < TempList[i].size(); k++)
                    {
                    m_angle_members.insert({TempList[i][j], h_tag.data[i], TempList[i][k]});
                    }
                }
            }
        TempList[i].clear();
        }

    //std::cout << "timestep=" << timestep;
    std::set<std::vector<unsigned int> > myList;
    for (const auto& vec : m_angle_members) 
        {
        assert(vec[0] <= m_pdata->getMaximumTag());
        assert(vec[1] <= m_pdata->getMaximumTag());
        assert(vec[2] <= m_pdata->getMaximumTag());

        // transform a, b, and c into indices into the particle data arrays
        // MEM TRANSFER: 6 ints
        unsigned int idx_a = h_rtag.data[vec[0]];
        unsigned int idx_b = h_rtag.data[vec[1]];
        unsigned int idx_c = h_rtag.data[vec[2]];

        // throw an error if this angle is incomplete
        if (idx_a == NOT_LOCAL || idx_b == NOT_LOCAL || idx_c == NOT_LOCAL)
            {
            this->m_exec_conf->msg->error()
                << "angle.harmonic: angle " << vec[0] << " " << vec[1] << " "
                << vec[2] << " incomplete." << endl
                << endl;
            throw std::runtime_error("Error in angle calculation");
            }

        assert(idx_a < m_pdata->getN() + m_pdata->getNGhosts());
        assert(idx_b < m_pdata->getN() + m_pdata->getNGhosts());
        assert(idx_c < m_pdata->getN() + m_pdata->getNGhosts());

        // calculate d\vec{r}
        Scalar3 dab;
        dab.x = h_pos.data[idx_a].x - h_pos.data[idx_b].x;
        dab.y = h_pos.data[idx_a].y - h_pos.data[idx_b].y;
        dab.z = h_pos.data[idx_a].z - h_pos.data[idx_b].z;

        Scalar3 dcb;
        dcb.x = h_pos.data[idx_c].x - h_pos.data[idx_b].x;
        dcb.y = h_pos.data[idx_c].y - h_pos.data[idx_b].y;
        dcb.z = h_pos.data[idx_c].z - h_pos.data[idx_b].z;

        Scalar3 dac;
        dac.x = h_pos.data[idx_a].x - h_pos.data[idx_c].x; // used for the 1-3 JL interaction
        dac.y = h_pos.data[idx_a].y - h_pos.data[idx_c].y;
        dac.z = h_pos.data[idx_a].z - h_pos.data[idx_c].z;

        // apply minimum image conventions to all 3 vectors
        dab = box.minImage(dab);
        dcb = box.minImage(dcb);
        dac = box.minImage(dac);

        // on paper, the formula turns out to be: F = K*\vec{r} * (r_0/r - 1)
        // FLOPS: 14 / MEM TRANSFER: 2 Scalars

        // FLOPS: 42 / MEM TRANSFER: 6 Scalars
        Scalar rsqab = dab.x * dab.x + dab.y * dab.y + dab.z * dab.z;
        Scalar rab = sqrt(rsqab);
        Scalar rsqcb = dcb.x * dcb.x + dcb.y * dcb.y + dcb.z * dcb.z;
        Scalar rcb = sqrt(rsqcb);

        Scalar hab = rab - 0.5*(h_diameter.data[idx_a]+h_diameter.data[idx_b]);
        Scalar hcb = rcb - 0.5*(h_diameter.data[idx_c]+h_diameter.data[idx_b]);
        if(hab>0.1 || hcb>0.1) myList.insert(vec);

        Scalar c_abbc = dab.x * dcb.x + dab.y * dcb.y + dab.z * dcb.z;
        c_abbc /= rab * rcb;
        //std::cout << "," << (int)(180.0/3.14*acos(c_abbc));
        if (c_abbc > 1.0)
            c_abbc = 1.0;
        if (c_abbc < -1.0)
            c_abbc = -1.0;

        Scalar s_abbc = sqrt(1.0 - c_abbc * c_abbc);
        if (s_abbc < SMALL)
            s_abbc = SMALL;
        s_abbc = 1.0 / s_abbc;

        // actually calculate the force
        unsigned int angle_type = m_angle_data->getTypeByIndex(0);
        Scalar dth = acos(c_abbc) - m_t_0[angle_type];
        Scalar tk = m_K[angle_type] * dth;
        Scalar a = -1.0 * tk * s_abbc;
        Scalar a11 = a * c_abbc / rsqab;
        Scalar a12 = -a / (rab * rcb);
        Scalar a22 = a * c_abbc / rsqcb;

        //std::cout << "timestep=" << timestep << "," << a11 << "," << a12 << "," << a22 << std::endl;
        Scalar fab[3], fcb[3];

        fab[0] = a11 * dab.x + a12 * dcb.x;
        fab[1] = a11 * dab.y + a12 * dcb.y;
        fab[2] = a11 * dab.z + a12 * dcb.z;

        fcb[0] = a22 * dcb.x + a12 * dab.x;
        fcb[1] = a22 * dcb.y + a12 * dab.y;
        fcb[2] = a22 * dcb.z + a12 * dab.z;

        // compute 1/3 of the energy, 1/3 for each atom in the angle
        Scalar angle_eng = (tk * dth) * Scalar(1.0 / 6.0);

        // compute 1/3 of the virial, 1/3 for each atom in the angle
        // upper triangular version of virial tensor
        Scalar angle_virial[6];
        angle_virial[0] = Scalar(1. / 3.) * (dab.x * fab[0] + dcb.x * fcb[0]);
        angle_virial[1] = Scalar(1. / 3.) * (dab.y * fab[0] + dcb.y * fcb[0]);
        angle_virial[2] = Scalar(1. / 3.) * (dab.z * fab[0] + dcb.z * fcb[0]);
        angle_virial[3] = Scalar(1. / 3.) * (dab.y * fab[1] + dcb.y * fcb[1]);
        angle_virial[4] = Scalar(1. / 3.) * (dab.z * fab[1] + dcb.z * fcb[1]);
        angle_virial[5] = Scalar(1. / 3.) * (dab.z * fab[2] + dcb.z * fcb[2]);

        // Now, apply the force to each individual atom a,b,c, and accumulate the energy/virial
        // do not update ghost particles
        if (idx_a < m_pdata->getN())
            {
            h_force.data[idx_a].x += fab[0];
            h_force.data[idx_a].y += fab[1];
            h_force.data[idx_a].z += fab[2];
            h_force.data[idx_a].w += angle_eng;
            for (int j = 0; j < 6; j++)
                h_virial.data[j * virial_pitch + idx_a] += angle_virial[j];
            }

        if (idx_b < m_pdata->getN())
            {
            h_force.data[idx_b].x -= fab[0] + fcb[0];
            h_force.data[idx_b].y -= fab[1] + fcb[1];
            h_force.data[idx_b].z -= fab[2] + fcb[2];
            h_force.data[idx_b].w += angle_eng;
            for (int j = 0; j < 6; j++)
                h_virial.data[j * virial_pitch + idx_b] += angle_virial[j];
            }

        if (idx_c < m_pdata->getN())
            {
            h_force.data[idx_c].x += fcb[0];
            h_force.data[idx_c].y += fcb[1];
            h_force.data[idx_c].z += fcb[2];
            h_force.data[idx_c].w += angle_eng;
            for (int j = 0; j < 6; j++)
                h_virial.data[j * virial_pitch + idx_c] += angle_virial[j];
            }
        }
    //std::cout << std::endl;

    for (auto& vec : myList) m_angle_members.erase(vec);
    myList.clear();

    //std::cout << "t=" << timestep << ",angles after addition=" << m_angle_members.size() << std::endl;
    }

namespace detail
    {
void export_ConstraintAngleForceCompute(pybind11::module& m)
    {
    pybind11::class_<ConstraintAngleForceCompute,
                     ForceCompute,
                     std::shared_ptr<ConstraintAngleForceCompute>>(m, "ConstraintAngleForceCompute")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>>())
        .def("setParams", &ConstraintAngleForceCompute::setParamsPython)
        .def("getParams", &ConstraintAngleForceCompute::getParams);
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
