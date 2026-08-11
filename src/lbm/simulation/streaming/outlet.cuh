#pragma once

#include "../../../config/config.h"
#include "../../../config/momentConfig.cuh"
#include "../../../config/simulationConfig.cuh"
#include "../../../config/stencilConfig.cuh"

#include "../../../core/indexing.cuh"

#include "../equations/f.cuh"

__device__ void outletEast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);