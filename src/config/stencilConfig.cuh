#pragma once

#include <math.h>

#include "config.h"
#include "gridConfig.cuh"

// Stencil definition
namespace D2Q9
{
    constexpr int D = 2;
    constexpr int Q = 9;

    constexpr int momNum = 6; // Moment number for D2Q9 is set at 6. Look at momentConfig.cuh for more details.

    constexpr int momSize = Geometry::NX * Geometry::NY;
    constexpr int momByteSize = Geometry::NX * Geometry::NY * momNum * sizeof(varType);

    const varType a_s = sqrtf(3);
    constexpr varType a_s2 = 3.f;
    constexpr varType a_s4 = 9.f;
    constexpr varType inv_as2 = 1.f / a_s2;

    __host__ __device__ __forceinline__ varType w(int i)
    {
        switch (i)
        {
        case 0:
            return 4.f / 9.f;
        case 1:
        case 2:
        case 3:
        case 4:
            return 1.f / 9.f;
        case 5:
        case 6:
        case 7:
        case 8:
            return 1.f / 36.f;
        default:
            return 0.f;
        }
    }

    __host__ __device__ __forceinline__ varType cx(int i)
    {
        switch (i)
        {
        case 0:
            return 0.f;
        case 1:
            return 1.f;
        case 2:
            return 0.f;
        case 3:
            return -1.f;
        case 4:
            return 0.f;
        case 5:
            return 1.f;
        case 6:
            return -1.f;
        case 7:
            return -1.f;
        case 8:
            return 1.f;
        default:
            return 0.f;
        }
    }

    __host__ __device__ __forceinline__ varType cy(int i)
    {
        switch (i)
        {
        case 0:
            return 0.f;
        case 1:
            return 0.f;
        case 2:
            return 1.f;
        case 3:
            return 0.f;
        case 4:
            return -1.f;
        case 5:
            return 1.f;
        case 6:
            return 1.f;
        case 7:
            return -1.f;
        case 8:
            return -1.f;
        default:
            return 0.f;
        }
    }

    __host__ __device__ __forceinline__ int income(int i)
    {
        switch (i)
        {
        case 0:
            return 0;
        case 1:
            return 3;
        case 2:
            return 4;
        case 3:
            return 1;
        case 4:
            return 2;
        case 5:
            return 7;
        case 6:
            return 8;
        case 7:
            return 5;
        default:
            return 6;
        }
    }
};