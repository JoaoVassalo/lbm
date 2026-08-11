#include "outlet.cuh"

__device__ void outletEast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
{
    maskType mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    // DIREÇÃO 0 (SEMPRE VERDADEIRA)

    int indexFrom = fromId(x, y, 0);

    varType fi = f(indexFrom, 0, mom_in);

    rhoI += fi;

    mxxI += fi * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += fi * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += fi * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);

    // OUTRAS DIREÇÕES
    while (mask)
    {
        int i = __ffs(mask);

        indexFrom = fromId(x, y, i);

        varType fi = f(indexFrom, i, mom_in);

        rhoI += fi;

        mxxI += fi * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += fi * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += fi * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    int neighbourId = fromId(x, y, (int)1);

    varType rho = mom_in[momIdx<momId::rho>(neighbourId)];
    varType ux = mom_in[momIdx<momId::ux>(neighbourId)];
    varType uy = mom_in[momIdx<momId::uy>(neighbourId)];

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(-rho - static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(3) * rho * ux) / (static_cast<varType>(6) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI + rho * uy) / (static_cast<varType>(3) * rho);
    mom_out[momIdx<momId::myy>(index)] = (static_cast<varType>(6) * myyI * rhoI) / (static_cast<varType>(5) * rho);
}