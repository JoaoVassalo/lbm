
#include "config/gridConfig.cuh"
#include "config/momentConfig.cuh"

#include "lbm/init/initialization.cuh"
#include "lbm/simulation/simulation.cuh"

int main()
{
    Grid2D grid;
    moments sim;

    initialization(sim, grid);

    simulation(sim, grid);
}