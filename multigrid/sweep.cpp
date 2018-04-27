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
#include "coarsen.h"

int main()
{
    DTMatlabDataFile inputFile("MyInput.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);
    DTMesh2D u0;
    Read(inputFile, "u", u0);

    DTMesh2DGrid grid = f.Grid();
    double h = grid.dx();
    std::vector<double> residuals;
    residuals.clear();

    DTMutableDoubleArray b = f.DoubleData().Copy();

    double scale = 1/(h*h);
    DTMutableDoubleArray r(b.m(),b.n());

    DTMatlabDataFile outputFile("OutputAll.mat", DTFile::NewReadWrite);
    DTMutableDoubleArray u = u0.DoubleData().Copy();

    sweep(u, b, 0.8, h, 100);
    sparseSolve();
    DTMatlabDataFile outputSweep("OutputSweep.mat",DTFile::NewReadWrite);

    DTMatlabDataFile referFile("ExactSolution.mat", DTFile::ReadOnly);
    DTMesh2D xexact;
    Read(referFile, "u", xexact);
    DTMutableDoubleArray x_exact = xexact.DoubleData().Copy();

    double err = 10000;
    for(int i = 1; i < x_exact.m()-1; i++)
    {
        for(int j = 1; j< x_exact.n()-1; ++j)
        {
            double err_norm = abs(x_exact(i,j) - u(i,j));
            if(err> err_norm)
            {
                err = err_norm;
            }
        }
    }
    std::cout << "eror norm:" << err << std::endl;
    Write(outputSweep, "u_final", u);
//    std::cout << "u"<< std::endl;
//    u.pall();
//    std::cout << "exact"<< std::endl;x_exact.pall();
//    std::cout << "b" << std::endl;
//    b.pall();

    double res = residual(u, b, scale, r);
//    r.pall();

    return 1;
}