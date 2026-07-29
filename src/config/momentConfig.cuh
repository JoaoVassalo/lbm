#pragma once

#include "config.h"

struct moments
{
    varType *momA;
    varType *momB;
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