//~ ########## Created by the Rheoinformatic research group ##########
//~ HOOMD-blue:
// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.
//~
//~ This file:
//~ Created by Rob Campbell (2022) [RHEOINF]

#include "AnisoPotentialPair.h"
#include "EvaluatorPairRotation.h"

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
inline void export_AnisoPotentialPair<EvaluatorPairRotation>(pybind11::module& m,
                                                           const std::string& name)
    {
    pybind11::class_<AnisoPotentialPair<EvaluatorPairRotation>,
                     ForceCompute,
                     std::shared_ptr<AnisoPotentialPair<EvaluatorPairRotation>>>
        anisopotentialpair(m, name.c_str());
    anisopotentialpair
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>>())
        .def("setParams", &AnisoPotentialPair<EvaluatorPairRotation>::setParamsPython)
        .def("getParams", &AnisoPotentialPair<EvaluatorPairRotation>::getParamsPython)
        .def("setMu", &AnisoPotentialPair<EvaluatorPairRotation>::setShapePython)
        .def("getMu", &AnisoPotentialPair<EvaluatorPairRotation>::getShapePython)
        .def("setRCut", &AnisoPotentialPair<EvaluatorPairRotation>::setRCutPython)
        .def("getRCut", &AnisoPotentialPair<EvaluatorPairRotation>::getRCut)
        .def_property("mode",
                      &AnisoPotentialPair<EvaluatorPairRotation>::getShiftMode,
                      &AnisoPotentialPair<EvaluatorPairRotation>::setShiftModePython)
        .def("getTypeShapesPy", &AnisoPotentialPair<EvaluatorPairRotation>::getTypeShapesPy);
    }

void export_AnisoPotentialPairRotation(pybind11::module& m)
    {
    export_AnisoPotentialPair<EvaluatorPairRotation>(m, "AnisoPotentialPairRotation");
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
