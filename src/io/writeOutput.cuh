#pragma once

#include "vtk.cuh"

#include "../config/config.h"
#include "../config/momentConfig.cuh"
#include "../config/gridConfig.cuh"

#include <string>

__host__ void writeOutput(varType *mom, varType *mom_host, Grid2D grid, size_t t, size_t momByteSize, const std::string path);