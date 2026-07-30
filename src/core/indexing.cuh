#pragma once

#include "../config/gridConfig.cuh"
#include "../config/stencilConfig.cuh"

struct coord
{
    int x;
    int y;
};

__device__ __forceinline__ coord coordinates()
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    return {x, y};
}

__device__ __forceinline__ int gridId()
{
    auto [x, y] = coordinates();
    return y * Geometry::NX + x;
}

__host__ __device__ __forceinline__ int gridId(int x, int y)
{
    return y * Geometry::NX + x;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int x, int y)
{
    return gridId(x, y) * D2Q9::momNum + I;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    return index * D2Q9::momNum + I;
}

__device__ __forceinline__ int fromId(int x, int y, int i)
{
    int x_from = x - (int)D2Q9::cx(i);
    int y_from = y - (int)D2Q9::cy(i);

    if (x_from < 0)
        x_from += Geometry::NX;
    if (x_from >= Geometry::NX)
        x_from -= Geometry::NX;

    if (y_from < 0)
        y_from += Geometry::NY;
    if (y_from >= Geometry::NY)
        y_from -= Geometry::NY;

    return gridId(x_from, y_from);
}