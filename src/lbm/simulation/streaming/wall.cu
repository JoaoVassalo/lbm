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

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI + static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega - static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega - static_cast<varType>(14) * ux - physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux - static_cast<varType>(14) * uy - physics::omega * uy - static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux - rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(static_cast<varType>(7) * rho + static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI + static_cast<varType>(18) * myyI * rhoI - static_cast<varType>(7) * rho * ux - static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    mom_out[momIdx<momId::myy>(index)] = -(static_cast<varType>(2) * (-rho + static_cast<varType>(6) * mxyI * rhoI - static_cast<varType>(9) * myyI * rhoI + rho * ux - static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
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

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI - static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega + static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega + static_cast<varType>(14) * ux + physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux - static_cast<varType>(14) * uy - physics::omega * uy + static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(static_cast<varType>(2) * (-rho - static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux + rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(7) * rho - static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI - static_cast<varType>(18) * myyI * rhoI - static_cast<varType>(7) * rho * ux + static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    mom_out[momIdx<momId::myy>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(9) * myyI * rhoI + rho * ux + static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
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

    varType ux = static_cast<varType>(physics::u_max);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI - static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega + static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega - static_cast<varType>(14) * ux - physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux + static_cast<varType>(14) * uy + physics::omega * uy + static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux + rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(7) * rho - static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI - static_cast<varType>(18) * myyI * rhoI + static_cast<varType>(7) * rho * ux - static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    mom_out[momIdx<momId::myy>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(9) * myyI * rhoI - rho * ux - static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
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

    varType ux = static_cast<varType>(physics::u_max);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(7) * mxyI * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI - static_cast<varType>(7) * mxyI * physics::omega * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega + static_cast<varType>(14) * ux + physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux + static_cast<varType>(14) * uy + physics::omega * uy - static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(static_cast<varType>(2) * (-rho - static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux - rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(static_cast<varType>(7) * rho + static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI + static_cast<varType>(18) * myyI * rhoI + static_cast<varType>(7) * rho * ux + static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    mom_out[momIdx<momId::myy>(index)] = -(static_cast<varType>(2) * (-rho + static_cast<varType>(6) * mxyI * rhoI - static_cast<varType>(9) * myyI * rhoI - rho * ux + static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
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

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega - static_cast<varType>(3) * uy - static_cast<varType>(3) * physics::omega * uy + static_cast<varType>(6) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (static_cast<varType>(6) * mxxI * rhoI) / (static_cast<varType>(5) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI - rho * ux) / (static_cast<varType>(3) * rho);
    mom_out[momIdx<momId::myy>(index)] = -(-rho - static_cast<varType>(9) * myyI * rhoI - static_cast<varType>(3) * rho * uy) / (static_cast<varType>(6) * rho);
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

    varType ux = static_cast<varType>(physics::u_max);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega + static_cast<varType>(3) * uy + static_cast<varType>(3) * physics::omega * uy + static_cast<varType>(6) * physics::omega * uy * uy);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = (static_cast<varType>(6) * mxxI * rhoI) / (static_cast<varType>(5) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI + rho * ux) / (static_cast<varType>(3) * rho);
    mom_out[momIdx<momId::myy>(index)] = -(-rho - static_cast<varType>(9) * myyI * rhoI + static_cast<varType>(3) * rho * uy) / (static_cast<varType>(6) * rho);
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

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega - static_cast<varType>(3) * ux - static_cast<varType>(3) * physics::omega * ux + static_cast<varType>(6) * physics::omega * ux * ux);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(-rho - static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(3) * rho * ux) /
                                         (static_cast<varType>(6) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI - rho * uy) /
                                         (static_cast<varType>(3) * rho);
    mom_out[momIdx<momId::myy>(index)] = (static_cast<varType>(6) * myyI * rhoI) /
                                         (static_cast<varType>(5) * rho);
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

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega + static_cast<varType>(3) * ux + static_cast<varType>(3) * physics::omega * ux + static_cast<varType>(6) * physics::omega * ux * ux);

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = -(-rho - static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(3) * rho * ux) / (static_cast<varType>(6) * rho);
    mom_out[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI + rho * uy) / (static_cast<varType>(3) * rho);
    mom_out[momIdx<momId::myy>(index)] = (static_cast<varType>(6) * myyI * rhoI) / (static_cast<varType>(5) * rho);
}