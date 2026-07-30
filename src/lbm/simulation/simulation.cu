#include "simulation.cuh"

__host__ void step(int t, varType *mom_in, varType *mom_out, Grid2D &grid)
{
    streamCollide<<<blockNumber, block>>>(mom_in, mom_out, grid);
}

__host__ void simulation(moments &sim, Grid2D &grid)
{
    auto t1 = std::chrono::high_resolution_clock::now;

    for (int t = 0; t < timeConfig::tf; t++)
    {
        if (t & 1)
            step(t, sim.momA, sim.momB, grid);
        else
            step(t, sim.momB, sim.momA, grid);
    }

    auto t2 = std::chrono::high_resolution_clock::now;
}