#include "initialization.cuh"

__host__ void initialization(moments &sim, Grid2D &grid)
{
    constexpr size_t momBytesize = D2Q9::momByteSize;

    cudaMalloc((void **)&sim.momA, momBytesize);
    cudaMalloc((void **)&sim.momB, momBytesize);

    initDomain<<<blockNumber, block>>>(sim.momA, sim.momB);

    cudaDeviceSynchronize();

    sim.mom_host = (varType *)malloc(momBytesize);

    cudaMemcpy(sim.mom_host, sim.momA, momBytesize, cudaMemcpyDeviceToHost);

    constexpr size_t maskBytesize = Grid2D::maskByteSize;
    constexpr size_t nodeBytesize = Grid2D::nodeByteSize;

    cudaMalloc((void **)&grid.mask, maskBytesize);
    cudaMalloc((void **)&grid.node, momBytesize);

    initGrid<<<blockNumber, block>>>(grid.mask, grid.node);
}