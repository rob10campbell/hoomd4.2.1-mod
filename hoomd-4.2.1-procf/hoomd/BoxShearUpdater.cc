//~ ########## Created by the PRO-CF research group ##########
//~ HOOMD-blue:
// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.
//~
//~ This file:
//~ Written by Dr. Deepak Mangal 
//~ Documentation by Rob Campbell (2024)

// ########## Created by PRO-CF //~ [PROCF2024] ##########

/*! \file BoxShearUpdater.cc
    \brief Defines the BoxShearUpdater class
*/

#include "BoxShearUpdater.h"

#include <iostream>
#include <math.h>
#include <stdexcept>

using namespace std;

namespace hoomd
    {
/*! \param sysdef System definition containing the particle data to set the box size on
    \param Lx length of the x dimension over time
    \param Ly length of the y dimension over time
    \param Lz length of the z dimension over time

    The default setting is to scale particle positions along with the box.
*/

BoxShearUpdater::BoxShearUpdater(std::shared_ptr<SystemDefinition> sysdef,
                                   std::shared_ptr<Trigger> trigger,
                                   std::shared_ptr<Variant> erate,
                                   Scalar deltaT,
                                   bool flip)
    : Updater(sysdef, trigger), m_erate(erate), m_deltaT(deltaT), m_flip(flip)
    {
    assert(m_pdata);
    assert(m_erate);
    m_exec_conf->msg->notice(5) << "Constructing BoxShearUpdater" << endl;
    }

BoxShearUpdater::~BoxShearUpdater()
    {
    m_exec_conf->msg->notice(5) << "Destroying BoxShearUpdater" << endl;
    }

void BoxShearUpdater::update(uint64_t timestep)
    {
    Updater::update(timestep);
    m_exec_conf->msg->notice(10) << "Box shear update" << endl;

    BoxDim cur_box = m_pdata->getGlobalBox();
    Scalar L_Y = cur_box.getL().y;

    Scalar cur_erate = (*m_erate)(timestep);
    Scalar3 new_L = cur_box.getL();
    Scalar xy = cur_box.getTiltFactorXY() + cur_erate * m_deltaT;
    if(m_flip)
        {
        if(xy>Scalar(0.50)) xy -= Scalar(1.0);
        else if(xy<Scalar(-0.50)) xy += Scalar(1.0);
        }
    Scalar xz = cur_box.getTiltFactorXZ();
    Scalar yz = cur_box.getTiltFactorYZ();
    BoxDim new_box = BoxDim(new_L);
    new_box.setTiltFactors(xy, xz, yz);
    m_pdata->setGlobalBox(new_box);

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<int3> h_image(m_pdata->getImages(),
                              access_location::host,
                              access_mode::readwrite);

    const BoxDim& local_box = m_pdata->getBox();
    for (unsigned int i = 0; i < m_pdata->getN(); i++)
        {
        int img0 = h_image.data[i].y; 
        local_box.wrap(h_pos.data[i], h_image.data[i]);
        img0 -= h_image.data[i].y;
        h_vel.data[i].x += (img0 * cur_erate * L_Y);
        }
    }

namespace detail
    {
void export_BoxShearUpdater(pybind11::module& m)
    {
    pybind11::class_<BoxShearUpdater, Updater, std::shared_ptr<BoxShearUpdater>>(
        m,
        "BoxShearUpdater")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<Trigger>,
                            std::shared_ptr<Variant>,Scalar, bool>())
        .def_property("erate", &BoxShearUpdater::getRate, &BoxShearUpdater::setRate)
        .def_property("deltaT", &BoxShearUpdater::getdeltaT, &BoxShearUpdater::setdeltaT)
        .def_property("flip", &BoxShearUpdater::getFlip, &BoxShearUpdater::setFlip);
    }

    } // end namespace detail

    } // end namespace hoomd
