#pragma once

#include "../config/config.h"
#include "../config/gridConfig.cuh"
#include "../config/momentConfig.cuh"

#include "../core/indexing.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>

__host__ void writeVTI(size_t t, const std::string &outDir, varType *mom_host);