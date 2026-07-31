#pragma once

#include "../../../config/config.h"
#include "../../../config/gridConfig.cuh"
#include "../../../config/simulationConfig.cuh"
#include "../../../config/stencilConfig.cuh"

#include "../../../core/indexing.cuh"

#include "../equations/f.cuh"

__device__ void wallSouthwest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallSoutheast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallNorthwest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallNortheast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallSouth(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallNorth(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallWest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);

__device__ void wallEast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index);