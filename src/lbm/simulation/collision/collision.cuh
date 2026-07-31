#pragma once

#include "../../../config/config.h"
#include "../../../config/momentConfig.cuh"
#include "../../../config/simulationConfig.cuh"

#include "../../../core/indexing.cuh"

__device__ void collide(varType *mom, int index);