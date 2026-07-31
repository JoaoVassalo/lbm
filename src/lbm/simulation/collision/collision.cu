#include "collision.cuh"

__device__ void collide(varType *mom, int index)
{
    const varType ux = mom[momIdx<momId::ux>(index)];
    const varType uy = mom[momIdx<momId::uy>(index)];

    const varType mxx = mom[momIdx<momId::mxx>(index)];
    const varType myy = mom[momIdx<momId::myy>(index)];
    const varType mxy = mom[momIdx<momId::mxy>(index)];

    mom[momIdx<momId::mxx>(index)] = ((1.0f - physics::omega) * mxx + physics::omega * ux * ux);
    mom[momIdx<momId::myy>(index)] = ((1.0f - physics::omega) * myy + physics::omega * uy * uy);
    mom[momIdx<momId::mxy>(index)] = ((1.0f - physics::omega) * mxy + physics::omega * ux * uy);
}