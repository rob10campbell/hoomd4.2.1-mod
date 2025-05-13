//~ ########## Created by the Rheoinformatic research group ##########
////~ HOOMD-blue:
//// Copyright (c) 2009-2022 The Regents of the University of Michigan.
//// Part of HOOMD-blue, released under the BSD 3-Clause License.
////~
////~ This file:
/////~ Written by Rob Campbell, Paniz Haghighi

#ifndef __STRETCH_MAP_H__
#define __STRETCH_MAP_H__

#pragma once

#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include "hoomd/GlobalArray.h"
#include "hoomd/HOOMDMath.h"
#include <iostream>
#include <memory>

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

// Declare the stretch_map as an extern variable
namespace hoomd
{
    extern std::map<std::pair<unsigned int, unsigned int>, Scalar3> stretch_map_temp;

    }
    #endif //__STRETCH_MAP_H__
