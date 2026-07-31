#pragma once

#include "gridConfig.cuh"

// Defining block size and number for kernel initialization, in this branch we use this part mainly for the initialization part.
#define BX 32
#define BY 16
#define GX (Geometry::NX / BX)
#define GY (Geometry::NY / BY)

inline dim3 block(BX, BY);
inline dim3 blockNumber(GX, GY);