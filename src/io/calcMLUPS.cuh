#pragma once

#include "../config/gridConfig.cuh"
#include "../config/simulationConfig.cuh"

#include <chrono>
#include <cstdio>

__host__ __forceinline__ void calcMLUPS(std::chrono::high_resolution_clock::time_point t1, std::chrono::high_resolution_clock::time_point t2)
{
    const auto elapsed = std::chrono::duration<double>(t2 - t1).count();

    const double mlups = static_cast<double>(Geometry::NX) * static_cast<double>(Geometry::NY) * static_cast<double>(timeConfig::tf) /
                         (elapsed * 1e6);

    printf("MLUPS: %.3f\n", mlups);
}