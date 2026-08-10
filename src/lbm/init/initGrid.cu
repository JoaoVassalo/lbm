#include "initGrid.cuh"

__global__ void initGrid(maskType *mask, nodeType *node)
{
    auto [x, y] = coordinates();

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    const int index = gridId(x, y);

    node[index] = toNodeType(defBC(x, y));

    maskType m = 0u;

    for (int i = 1; i < D2Q9::Q; i++)
    {
        int xn = x - D2Q9::cx(i);
        int yn = y - D2Q9::cy(i);

        if (xn < 0 || xn >= Geometry::NX ||
            yn < 0 || yn >= Geometry::NY)
        {
            continue;
        }

        m |= (1u << (i - 1));
    }

    mask[index] = m;
}