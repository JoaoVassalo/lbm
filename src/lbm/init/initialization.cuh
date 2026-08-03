#pragma once

#include "../../config/config.h"
#include "../../config/CUDAConfig.cuh"
#include "../../config/gridConfig.cuh"
#include "../../config/momentConfig.cuh"
#include "../../config/outputConfig.h"
#include "../../config/stencilConfig.cuh"

#include "../../core/indexing.cuh"

#include "../../io/writeOutput.cuh"

#include "initDomain.cuh"
#include "initGrid.cuh"

__host__ void initialization(moments &sim, Grid2D &grid);