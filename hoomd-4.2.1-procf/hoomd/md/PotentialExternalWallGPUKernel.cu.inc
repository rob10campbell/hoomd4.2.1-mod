// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// See md/CMakeLists.txt for the source of these variables to be processed by CMake's
// configure_file().

// clang-format off
#include "hoomd/md/PotentialExternalGPU.cuh"
#include "hoomd/md/EvaluatorPair@_evaluator@.h"
#include "hoomd/md/EvaluatorWalls.h"

#define EVALUATOR_CLASS EvaluatorWalls<EvaluatorPair@_evaluator@>
// clang-format on

namespace hoomd
    {
namespace md
    {
namespace kernel
    {
template __attribute__((visibility("default"))) hipError_t
gpu_compute_potential_external_forces<EVALUATOR_CLASS>(
    const external_potential_args_t& external_potential_args,
    const typename EVALUATOR_CLASS::param_type* d_params,
    const typename EVALUATOR_CLASS::field_type* d_field);
    } // end namespace kernel
    } // end namespace md
    } // end namespace hoomd