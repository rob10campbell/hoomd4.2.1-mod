//~ ########## Created by the Rheoinformatic research group ##########
//~ HOOMD-blue:
// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.
//~
//~ This file:
//~ Created by Rob Campbell (2022) [RHEOINF]

#include "AnisoPotentialPair.h"
#include "EvaluatorPairMorseFrix.h"

namespace hoomd
    {
namespace md
    {
namespace detail
    {

// Template specification modeled after Dipole anisotropic pair potential. A specific
// template instance is needed since we expose the shape as just mu in Python
// when the default behavior exposes setting and getting the shape through
// 'shape'.
template<>
inline void export_AnisoPotentialPair<EvaluatorPairMorseFrix>(pybind11::module& m,
                                                           const std::string& name)
    {
    pybind11::class_<AnisoPotentialPair<EvaluatorPairMorseFrix>,
                     ForceCompute,
                     std::shared_ptr<AnisoPotentialPair<EvaluatorPairMorseFrix>>>
        anisopotentialpair(m, name.c_str());
    anisopotentialpair
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, //>() //~ add angular rigidity params [RHEOINF]
                            Scalar, 
                            Scalar, 
                            Scalar>()) //~ add angular rigidity params [RHEOINF]
        .def("setParams", &AnisoPotentialPair<EvaluatorPairMorseFrix>::setParamsPython)
        .def("getParams", &AnisoPotentialPair<EvaluatorPairMorseFrix>::getParamsPython)
        .def("setMu", &AnisoPotentialPair<EvaluatorPairMorseFrix>::setShapePython)
        .def("getMu", &AnisoPotentialPair<EvaluatorPairMorseFrix>::getShapePython)
        .def("setRCut", &AnisoPotentialPair<EvaluatorPairMorseFrix>::setRCutPython)
        .def("getRCut", &AnisoPotentialPair<EvaluatorPairMorseFrix>::getRCut)
        .def_property("mode",
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::getShiftMode,
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::setShiftModePython)
        .def_property("K", 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::getK, 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::setK) //~ add K [RHEOINF]
        .def_property("w", 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::getW, 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::setW) //~ add w [RHEOINF]
        .def_property("theta_bar", 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::getTheta, 
                      &AnisoPotentialPair<EvaluatorPairMorseFrix>::setTheta) //~ add theta_bar [RHEOINF]
        .def("getTypeShapesPy", &AnisoPotentialPair<EvaluatorPairMorseFrix>::getTypeShapesPy);
    }

void export_AnisoPotentialPairMorseFrix(pybind11::module& m)
    {
    export_AnisoPotentialPair<EvaluatorPairMorseFrix>(m, "AnisoPotentialPairMorseFrix");
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
