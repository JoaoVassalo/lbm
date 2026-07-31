#pragma once

#include <numbers>

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

    constexpr varType a_s = static_cast<varType>(std::numbers::sqrt3_v<double>);
    constexpr varType a_s2 = static_cast<varType>(3);
    constexpr varType a_s4 = static_cast<varType>(9);
    constexpr varType inv_as2 = static_cast<varType>(1.0 / a_s2);

    __host__ __device__ __forceinline__ float w(int i)
    {
        switch (i)
        {
        case 0:
            return static_cast<varType>(4.0 / 9.0);
        case 1:
        case 2:
        case 3:
        case 4:
            return static_cast<varType>(1.0 / 9.0);
        case 5:
        case 6:
        case 7:
        case 8:
            return static_cast<varType>(1.0 / 36.0);
        default:
            return static_cast<varType>(0.0);
        }
    }

    __host__ __device__ __forceinline__ varType cx(int i)
    {
        switch (i)
        {
        case 0:
            return static_cast<varType>(0);
        case 1:
            return static_cast<varType>(1);
        case 2:
            return static_cast<varType>(0);
        case 3:
            return static_cast<varType>(-1);
        case 4:
            return static_cast<varType>(0);
        case 5:
            return static_cast<varType>(1);
        case 6:
            return static_cast<varType>(-1);
        case 7:
            return static_cast<varType>(-1);
        case 8:
            return static_cast<varType>(1);
        default:
            return static_cast<varType>(0);
        }
    }

    __host__ __device__ __forceinline__ float cy(int i)
    {
        switch (i)
        {
        case 0:
            return static_cast<varType>(0);
        case 1:
            return static_cast<varType>(0);
        case 2:
            return static_cast<varType>(1);
        case 3:
            return static_cast<varType>(0);
        case 4:
            return -static_cast<varType>(1);
        case 5:
            return static_cast<varType>(1);
        case 6:
            return static_cast<varType>(1);
        case 7:
            return -static_cast<varType>(1);
        case 8:
            return -static_cast<varType>(1);
        default:
            return static_cast<varType>(0);
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