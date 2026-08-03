#include "initDomain.cuh"

__global__ void initDomain(varType *momA, varType *momB)
{
    auto [x, y] = coordinates();

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    momA[momIdx<momId::rho>(x, y)] = 1.f;
    momA[momIdx<momId::ux>(x, y)] = 0.f;
    momA[momIdx<momId::uy>(x, y)] = 0.f;
    momA[momIdx<momId::mxx>(x, y)] = 0.f;
    momA[momIdx<momId::mxy>(x, y)] = 0.f;
    momA[momIdx<momId::myy>(x, y)] = 0.f;

    momB[momIdx<momId::rho>(x, y)] = 1.f;
    momB[momIdx<momId::ux>(x, y)] = 0.f;
    momB[momIdx<momId::uy>(x, y)] = 0.f;
    momB[momIdx<momId::mxx>(x, y)] = 0.f;
    momB[momIdx<momId::mxy>(x, y)] = 0.f;
    momB[momIdx<momId::myy>(x, y)] = 0.f;
}