#pragma once

#include "config.h"
#include "gridConfig.cuh"

struct moments
{
    varType *mom;
    varType *mom_host;
};

enum momId
{
    rho = 0,
    ux,
    uy,
    mxx,
    mxy,
    myy,
    count
};