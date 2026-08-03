#include "fluid.cuh"

__device__ void fluid(varType *mom_in, varType *mom_out, int x, int y, int index)
{
    varType rho = 0.f;
    varType ux = 0.f;
    varType uy = 0.f;
    varType mxx = 0.f;
    varType mxy = 0.f;
    varType myy = 0.f;

#pragma unroll
    for (int i = 0; i < D2Q9::Q; i++)
    {
        int indexFrom = fromId(x, y, i);
        varType fi = f(indexFrom, i, mom_in);

        rho += fi;

        ux += fi * D2Q9::cx(i);
        uy += fi * D2Q9::cy(i);

        mxx += fi * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxy += fi * (D2Q9::cx(i) * D2Q9::cy(i));
        myy += fi * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);
    }

    ux /= rho;
    uy /= rho;
    mxx /= rho;
    mxy /= rho;
    myy /= rho;

    mom_out[momIdx<momId::rho>(index)] = rho;
    mom_out[momIdx<momId::ux>(index)] = ux;
    mom_out[momIdx<momId::uy>(index)] = uy;
    mom_out[momIdx<momId::mxx>(index)] = mxx;
    mom_out[momIdx<momId::mxy>(index)] = mxy;
    mom_out[momIdx<momId::myy>(index)] = myy;
}