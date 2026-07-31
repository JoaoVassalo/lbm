#include "wall.cuh"

__device__ void wallSouthwest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (12.f * (-3.f - 3.f * mxxI + 7.f * mxyI - 3.f * myyI + 3.f * mxxI * physics::omega - 7.f * mxyI * physics::omega + 3.f * myyI * physics::omega) * rhoI) /
                  (-16.f - 9.f * physics::omega - 14.f * ux - physics::omega * ux + 15.f * physics::omega * ux * ux - 14.f * uy - physics::omega * uy - 9.f * physics::omega * ux * uy + 15.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI - 6.f * mxyI * rhoI + 2.f * rho * ux - rho * uy)) /
                                         (9.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI + 18.f * myyI * rhoI - 7.f * rho * ux - 7.f * rho * uy) /
                                         (27.f * rho);
    mom_out[momIdx<momId::myy>(index)] = -(2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI + rho * ux - 2.f * rho * uy)) /
                                         (9.f * rho);
}

__device__ void wallSoutheast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (12.f * (-3.f - 3.f * mxxI - 7.f * mxyI - 3.f * myyI + 3.f * mxxI * physics::omega + 7.f * mxyI * physics::omega + 3.f * myyI * physics::omega) * rhoI) /
                  (-16.f - 9.f * physics::omega + 14.f * ux + physics::omega * ux + 15.f * physics::omega * ux * ux - 14.f * uy - physics::omega * uy + 9.f * physics::omega * ux * uy + 15.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(2.f * (-rho - 9.f * mxxI * rhoI - 6.f * mxyI * rhoI + 2.f * rho * ux + rho * uy)) /
                                         (9.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI - 18.f * myyI * rhoI - 7.f * rho * ux + 7.f * rho * uy) /
                                         (27.f * rho);
    mom_out[momIdx<momId::myy>(index)] = (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI + rho * ux + 2.f * rho * uy)) /
                                         (9.f * rho);
}

__device__ void wallNorthwest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (12.f * (-3.f - 3.f * mxxI - 7.f * mxyI - 3.f * myyI + 3.f * mxxI * physics::omega + 7.f * mxyI * physics::omega + 3.f * myyI * physics::omega) * rhoI) /
                  (-16.f - 9.f * physics::omega - 14.f * ux - physics::omega * ux + 15.f * physics::omega * ux * ux + 14.f * uy + physics::omega * uy + 9.f * physics::omega * ux * uy + 15.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI + 6.f * mxyI * rhoI + 2.f * rho * ux + rho * uy)) /
                                         (9.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI - 18.f * myyI * rhoI + 7.f * rho * ux - 7.f * rho * uy) /
                                         (27.f * rho);
    mom_out[momIdx<momId::myy>(index)] = (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI - rho * ux - 2.f * rho * uy)) /
                                         (9.f * rho);
}

__device__ void wallNortheast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (12.f * (-3.f * rhoI - 3.f * mxxI * rhoI + 7.f * mxyI * rhoI - 3.f * myyI * rhoI + 3.f * mxxI * physics::omega * rhoI - 7.f * mxyI * physics::omega * rhoI + 3.f * myyI * physics::omega * rhoI)) /
                  (-16.f - 9.f * physics::omega + 14.f * ux + physics::omega * ux + 15.f * physics::omega * ux * ux + 14.f * uy + physics::omega * uy - 9.f * physics::omega * ux * uy + 15.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(2.f * (-rho - 9.f * mxxI * rhoI + 6.f * mxyI * rhoI + 2.f * rho * ux - rho * uy)) /
                                         (9.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI + 18.f * myyI * rhoI + 7.f * rho * ux + 7.f * rho * uy) /
                                         (27.f * rho);
    mom_out[momIdx<momId::myy>(index)] = -(2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI - rho * ux + 2.f * rho * uy)) /
                                         (9.f * rho);
}

__device__ void wallSouth(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (3.f * (-4.f * rhoI - 3.f * myyI * rhoI + 3.f * myyI * physics::omega * rhoI)) /
                  (-9.f - physics::omega - 3.f * uy - 3.f * physics::omega * uy + 6.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-6.f * mxyI * rhoI - rho * ux) / (3.f * rho);
    mom_out[momIdx<momId::myy>(index)] = -(-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
}

__device__ void wallNorth(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (3.f * (-4.f * rhoI - 3.f * myyI * rhoI + 3.f * myyI * physics::omega * rhoI)) /
                  (-9.f - physics::omega + 3.f * uy + 3.f * physics::omega * uy + 6.f * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-6.f * mxyI + rho * ux) / (3.f * rho);
    mom_out[momIdx<momId::myy>(index)] = -(-rho - 9.f * myyI * rhoI + 3.f * rho * uy) / (6.f * rho);
}

__device__ void wallWest(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (3.f * (-4.f * rhoI - 3.f * mxxI * rhoI + 3.f * mxxI * physics::omega * rhoI)) /
                  (-9.f - physics::omega - 3.f * ux - 3.f * physics::omega * ux + 6.f * physics::omega * ux * ux);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(-rho - 9.f * mxxI * rhoI - 3.f * rho * ux) /
                                         (6.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-6.f * mxyI * rhoI - rho * uy) /
                                         (3.f * rho);
    mom_out[momIdx<momId::myy>(index)] = (6.f * myyI * rhoI) /
                                         (5.f * rho);
}

__device__ void wallEast(int x, int y, varType *mom_in, varType *mom_out, Grid2D grid, int index)
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
        int bit = __ffs(mask);

        int i = bit + 1;

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

    varType ux = 0.f;
    varType uy = 0.f;

    varType rho = (3.f * (-4.f * rhoI - 3.f * mxxI * rhoI + 3.f * mxxI * physics::omega * rhoI)) /
                  (-9.f - physics::omega + 3.f * ux + 3.f * physics::omega * ux + 6.f * physics::omega * ux * ux);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(-rho - 9.f * mxxI * rhoI + 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-6.f * mxyI * rhoI + rho * uy) / (3.f * rho);
    mom_out[momIdx<momId::myy>(index)] = (6.f * myyI * rhoI) / (5.f * rho);
}