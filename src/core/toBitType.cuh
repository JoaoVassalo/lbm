#pragma once

#include "../config/config.h"

#include <cstdint>

__device__ __forceinline__ maskType toMaskType(auto num)
{
    return static_cast<maskType>(num);
}

__device__ __forceinline__ nodeType toNodeType(auto num)
{
    return static_cast<nodeType>(num);
}