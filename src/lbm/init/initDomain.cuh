#pragma once

#include "../../config/config.h"
#include "../../config/momentConfig.cuh"
#include "../../config/stencilConfig.cuh"

#include "../../core/indexing.cuh"

__global__ void initDomain(varType *momA, varType *momB);