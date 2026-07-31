#include "wall.cuh"

__device__ void wallSouthwest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid)
{
    int index = gridId(x, y);

    maskType mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    while (mask)
    {
        /* code */
    }
}

__device__ void wallSoutheast()
{
}

__device__ void wallNorthwest()
{
}

__device__ void wallNortheast()
{
}

__device__ void wallSouth()
{
}

__device__ void wallNorth()
{
}

__device__ void wallWest()
{
}

__device__ void wallEast()
{
}