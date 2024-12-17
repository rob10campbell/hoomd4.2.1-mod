// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "AnisoPotentialPair.h"
#include "EvaluatorPairDipole.h"

namespace hoomd
    {
namespace md
    {
namespace detail
    {

// Template specification for Dipole anisotropic pair potential. A specific
// template instance is needed since we expose the shape as just mu in Python
// when the default behavior exposes setting and getting the shape through
// 'shape'.
template<>
inline void export_AnisoPotentialPair<EvaluatorPairDipole>(pybind11::module& m,
                                                           const std::string& name)
    {
    pybind11::class_<AnisoPotentialPair<EvaluatorPairDipole>,
                     ForceCompute,
                     std::shared_ptr<AnisoPotentialPair<EvaluatorPairDipole>>>
        anisopotentialpair(m, name.c_str());
    anisopotentialpair
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<NeighborList>, //>() //~ add angular rigidity params [RHEOINF]
                            Scalar,
                            Scalar,
                            Scalar>()) //~ add angular rigidity params [RHEOINF]
        .def("setParams", &AnisoPotentialPair<EvaluatorPairDipole>::setParamsPython)
        .def("getParams", &AnisoPotentialPair<EvaluatorPairDipole>::getParamsPython)
        .def("setMu", &AnisoPotentialPair<EvaluatorPairDipole>::setShapePython)
        .def("getMu", &AnisoPotentialPair<EvaluatorPairDipole>::getShapePython)
        .def("setRCut", &AnisoPotentialPair<EvaluatorPairDipole>::setRCutPython)
        .def("getRCut", &AnisoPotentialPair<EvaluatorPairDipole>::getRCut)
        .def_property("mode",
                      &AnisoPotentialPair<EvaluatorPairDipole>::getShiftMode,
                      &AnisoPotentialPair<EvaluatorPairDipole>::setShiftModePython)
        .def_property("K",
                      &AnisoPotentialPair<EvaluatorPairDipole>::getK,
                      &AnisoPotentialPair<EvaluatorPairDipole>::setK) //~ add K [RHEOINF]
        .def_property("w",
                      &AnisoPotentialPair<EvaluatorPairDipole>::getW,
                      &AnisoPotentialPair<EvaluatorPairDipole>::setW) //~ add w [RHEOINF]
        .def_property("theta_bar",
                      &AnisoPotentialPair<EvaluatorPairDipole>::getTheta,
                      &AnisoPotentialPair<EvaluatorPairDipole>::setTheta) //~ add theta_bar [RHEOINF]
        .def("getTypeShapesPy", &AnisoPotentialPair<EvaluatorPairDipole>::getTypeShapesPy);
    }

void export_AnisoPotentialPairDipole(pybind11::module& m)
    {
    export_AnisoPotentialPair<EvaluatorPairDipole>(m, "AnisoPotentialPairDipole");
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
