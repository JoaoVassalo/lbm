#pragma once

#include "../../../config/config.h"
#include "../../../config/momentConfig.cuh"

#include "../../../core/indexing.cuh"

__device__ __forceinline__ varType f(int index, int i, varType *mom)
{
    return mom[momIdx<momId::rho>(index)] * D2Q9::w(i) *
           (static_cast<varType>(1.0) +
            D2Q9::a_s2 * mom[momIdx<momId::ux>(index)] * D2Q9::cx(i) +
            D2Q9::a_s2 * mom[momIdx<momId::uy>(index)] * D2Q9::cy(i) +
            D2Q9::a_s4 * static_cast<varType>(0.5) * mom[momIdx<momId::mxx>(index)] * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2) +
            D2Q9::a_s4 * mom[momIdx<momId::mxy>(index)] * (D2Q9::cx(i) * D2Q9::cy(i)) +
            D2Q9::a_s4 * static_cast<varType>(0.5) * mom[momIdx<momId::myy>(index)] * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2));
}