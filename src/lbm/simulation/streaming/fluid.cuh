#pragma once

#include "../../../config/config.h"
#include "../../../config/stencilConfig.cuh"

#include "../../../core/indexing.cuh"

#include "../equations/f.cuh"

__device__ void fluid(varType *mom_in, varType *mom_out, int x, int y, int index);