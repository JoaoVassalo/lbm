#pragma once

#include "../../config/config.h"
#include "../../config/gridConfig.cuh"
#include "../../config/stencilConfig.cuh"

#include "../../core/indexing.cuh"
#include "../../core/toBitType.cuh"

__global__ void initGrid(maskType *mask, nodeType *node);