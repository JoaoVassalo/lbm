#pragma once

#include "config.h"
#include "gridConfig.cuh"
#include "stencilConfig.cuh"

namespace timeConfig
{
    inline constexpr size_t tf = 1000;
    inline constexpr size_t printNumber = 100;
    inline constexpr size_t tInterval = tf / printNumber;
}

namespace physics
{
    inline constexpr int Re = 1000;
    inline constexpr varType u_max = static_cast<varType>(0.0256);
    inline constexpr varType delta_t = static_cast<varType>(1.0);
    inline constexpr varType ni = u_max * static_cast<varType>(Geometry::NY) / static_cast<varType>(Re);
    inline constexpr varType tau = ni * D2Q9::a_s2 + static_cast<varType>(0.5);
    inline constexpr varType omega = static_cast<varType>(1.0) / tau;
}