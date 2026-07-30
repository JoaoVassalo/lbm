#include "wall.cuh"

__device__ void wallSouthwest(int x, int y, varType *mom_in, varType *mom_out)
{
    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    for (int i = 0; i < D2Q9::Q; i++)
    {
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