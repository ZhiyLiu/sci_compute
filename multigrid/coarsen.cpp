#include "coarsen.h"
#include "stdlib.h"
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "DTMatlabDataFile.h"
#include "DTArguments.h"
#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTMesh2D.h"
#include <sstream>
#include "DTFunction2D.h"
#include "DTMesh2DGrid.h"
#include "DTSeriesMesh2D.h"
#include "time.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "sweep.h"
#include "sparseSolve.h"
#include "residual.h"

/********************************
 //    Test Coarsen & Refine
 *************************/
int main(){
    int dim = 9; // 9*9 include bdry should coarsen to 5*5
    double h = 0.025;
    DTPoint2D origin(0.0, 0.0);

    DTMutableDoubleArray data(dim, dim);
    for(int i = 0; i < dim; ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            double x = i* h;
            double y = j * h;
            data(i,j) = x + 5 * y;
        }
    }
    DTMesh2D fineGrid(DTMesh2DGrid(origin, h, h, dim, dim), data);

    DTMutableMesh2D coarseGrid(DTMesh2DGrid(origin, 2*h, 2*h,1+ dim/2, 1 + dim/2), DTMutableDoubleArray((int)(1+dim/2), (int)(1+dim/2)));
    coarsen(fineGrid, coarseGrid);

    DTMutableMesh2D refineGrid(DTMesh2DGrid(origin, h, h, dim, dim), DTMutableDoubleArray(dim,dim));
    refine(coarseGrid, refineGrid);
    DTMatlabDataFile outputCoarsen("Coarsen.mat", DTFile::NewReadWrite);
    Write(outputCoarsen, "fine", fineGrid);
    Write(outputCoarsen, "coarse", coarseGrid);
    Write(outputCoarsen, "refine", refineGrid);
    return 0;
}