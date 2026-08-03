#include "simulation.cuh"

__host__ void step(size_t t, varType *mom_in, varType *mom_out, varType *mom_host, Grid2D &grid)
{
    streamCollide<<<blockNumber, block>>>(mom_in, mom_out, grid);

    cudaDeviceSynchronize();
    if (t % timeConfig::tInterval == 0)
        writeOutput(mom_out, mom_host, grid, t, D2Q9::momByteSize, output::vtkPath);
}

__host__ void simulation(moments &sim, Grid2D &grid)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t t = 1; t < timeConfig::tf; t++)
    {
        if (t & 1)
            step(t, sim.momA, sim.momB, sim.mom_host, grid);
        else
            step(t, sim.momB, sim.momA, sim.mom_host, grid);
    }

    cudaDeviceSynchronize();

    auto t2 = std::chrono::high_resolution_clock::now();

    calcMLUPS(t1, t2);
}