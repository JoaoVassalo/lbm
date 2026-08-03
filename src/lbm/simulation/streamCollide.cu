#include "streamCollide.cuh"

__device__ void applyBoundary(varType *mom_in, varType *mom_out, Grid2D grid, int x, int y, int index)
{
    if (x == 0 && y == 0)
    { // Sudoeste
        wallSouthwest(x, y, mom_in, mom_out, grid, index);
    }
    else if (x == Geometry::NX - 1 && y == 0)
    { // Sudeste
        wallSoutheast(x, y, mom_in, mom_out, grid, index);
    }
    else if (x == 0 && y == Geometry::NY - 1)
    { // Noroeste
        wallNorthwest(x, y, mom_in, mom_out, grid, index);
    }
    else if (x == Geometry::NX - 1 && y == Geometry::NY - 1)
    { // Nordeste
        wallNortheast(x, y, mom_in, mom_out, grid, index);
    }
    else if (y == 0)
    { // Sul
        wallSouth(x, y, mom_in, mom_out, grid, index);
    }
    else if (y == Geometry::NY - 1)
    { // Norte
        wallNorth(x, y, mom_in, mom_out, grid, index);
    }
    else if (x == 0)
    { // Oeste
        wallWest(x, y, mom_in, mom_out, grid, index);
    }
    else if (x == Geometry::NX - 1)
    { // Leste
        wallEast(x, y, mom_in, mom_out, grid, index);
    }
    else
    { // Centro
        fluid(mom_in, mom_out, x, y, index);
    }
}

__global__ void streamCollide(varType *mom_in, varType *mom_out, Grid2D grid)
{
    auto [x, y] = coordinates();

    if (x >= Geometry::NX || y >= Geometry::NY || x < 0 || y < 0)
        return;

    int index = gridId();

    applyBoundary(mom_in, mom_out, grid, x, y, index);
    collide(mom_out, index);
}