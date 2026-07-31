#include "writeOutput.cuh"

__host__ void writeOutput(varType *mom, varType *mom_host, Grid2D grid, size_t t, size_t momByteSize, const std::string path)
{
    cudaDeviceSynchronize();
    cudaMemcpy(mom_host, mom, momByteSize, cudaMemcpyDeviceToHost);
    writeVTI(t, path, mom_host);

    std::cout << "Iteration " << t << std::endl;
    std::cout << std::endl;
}