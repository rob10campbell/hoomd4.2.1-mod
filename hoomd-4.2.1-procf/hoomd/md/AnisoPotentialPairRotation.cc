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
template void export_AnisoPotentialPair<EvaluatorPairRotation>(pybind11::module& m,
                                                         const std::string& name);

void export_AnisoPotentialPairRotation(pybind11::module& m)
    {
    export_AnisoPotentialPair<EvaluatorPairRotation>(m, "AnisoPotentialPairRotation");
    }
    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
