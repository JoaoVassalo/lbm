#include "streamCollide.cuh"

__device__ void applyBoundary(varType *mom_in, varType *mom_out, Grid2D grid, int x, int y)
{
    if (x == 0 && y == 0)
    { // Sudoeste
    }
    else if (x == Geometry::NX - 1 && y == 0)
    { // Sudeste
    }
    else if (x == 0 && y == Geometry::NY - 1)
    { // Noroeste
    }
    else if (x == Geometry::NX - 1 && y == Geometry::NY - 1)
    { // Nordeste
    }
    else if (y == 0)
    { // Sul
    }
    else if (y == Geometry::NY - 1)
    { // Norte
    }
    else if (x == 0)
    { // Oeste
    }
    else if (x == Geometry::NX - 1)
    { // Leste
    }
    else
    { // Centro
        fluid(mom_in, mom_out, x, y);
    }
}

__global__ void streamCollide(varType *mom_in, varType *mom_out, Grid2D grid)
{
    auto [x, y] = coordinates();

    if (x >= Geometry::NX || y >= Geometry::NY || x < 0 || y < 0)
        return;

    applyBoundary(mom_in, mom_out, grid, x, y);
    collide();
}