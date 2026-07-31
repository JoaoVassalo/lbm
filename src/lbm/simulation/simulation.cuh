#pragma once

#include <chrono>

#include "../../config/config.h"
#include "../../config/CUDAConfig.cuh"
#include "../../config/gridConfig.cuh"
#include "../../config/outputConfig.h"
#include "../../config/simulationConfig.cuh"
#include "../../config/momentConfig.cuh"

#include "../../io/calcMLUPS.cuh"
#include "../../io/writeOutput.cuh"

#include "streamCollide.cuh"

__host__ void simulation(moments &sim, Grid2D &grid);