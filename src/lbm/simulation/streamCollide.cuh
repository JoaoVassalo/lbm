#pragma once

#include "../../config/config.h"
#include "../../config/gridConfig.cuh"
#include "../../config/momentConfig.cuh"

#include "../../core/indexing.cuh"

#include "collision/collision.cuh"

#include "streaming/fluid.cuh"
#include "streaming/inlet.cuh"
#include "streaming/outlet.cuh"
#include "streaming/wall.cuh"

__global__ void streamCollide(varType *mom_in, varType *mom_out, Grid2D grid);