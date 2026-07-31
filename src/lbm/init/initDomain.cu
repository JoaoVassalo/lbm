#include "initDomain.cuh"

__global__ void initDomain(varType *momA, varType *momB)
{
    auto [x, y] = coordinates();

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    momB[momIdx<momId::rho>(x, y)] = 1.f;
    momB[momIdx<momId::ux>(x, y)] = 0.f;
    if (x == 0)
        momB[momIdx<momId::ux>(x, y)] = physics::u_max;
    momB[momIdx<momId::uy>(x, y)] = 0.f;
    momB[momIdx<momId::mxx>(x, y)] = 0.f;
    momB[momIdx<momId::mxy>(x, y)] = 0.f;
    momB[momIdx<momId::myy>(x, y)] = 0.f;
}