#pragma once

#include "config.h"

#include <cstdint>

namespace Geometry
{
    // Geometry definition
    constexpr int NX = 512;
    constexpr int NY = 128;
    constexpr int NZ = 0;
};

struct Grid2D
{
    maskType *mask;
    nodeType *node;

    static constexpr int maskByteSize = Geometry::NX * Geometry::NY * sizeof(maskType);
    static constexpr int nodeByteSize = Geometry::NX * Geometry::NY * sizeof(nodeType);
};

enum Boundary
{
    Center = 0,
    East,
    North,
    West,
    South,
    Northeast,
    Northwest,
    Southwest,
    Southeast
};

__host__ __device__ __forceinline__ Boundary defBC(int x, int y)
{
    if (x == 0 && y == 0)
    {
        return Boundary::Southwest;
    }
    else if (x == Geometry::NX - 1 && y == 0)
    {
        return Boundary::Southeast;
    }
    else if (x == 0 && y == Geometry::NY - 1)
    {
        return Boundary::Northwest;
    }
    else if (x == Geometry::NX - 1 && y == Geometry::NY - 1)
    {
        return Boundary::Northeast;
    }
    else if (y == 0)
    {
        return Boundary::South;
    }
    else if (y == Geometry::NY - 1)
    {
        return Boundary::North;
    }
    else if (x == 0)
    {
        return Boundary::West;
    }
    else if (x == Geometry::NX - 1)
    {
        return Boundary::East;
    }
    else
    {
        return Boundary::Center;
    }
}