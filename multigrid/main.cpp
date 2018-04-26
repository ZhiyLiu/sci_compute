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

void recurse(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double omega, double h, int Nb, int Na, DTPoint2D origin)
{
    // terminate criteria
//    if(h >= h_max)
    if(b.m() < 6 || b.n() < 6)
    {
        //solve it if the size excluding boundary is 3*3 or 2*2
//        sweep(u, b, omega, h, 1000);
        directSolve(b.m()-2, b.n()-2, h, b, u);
        return;
    }

    /****************************
     //    Sweep before
     *************************/
    sweep(u, b, omega, h, Nb);

    /****************************
     //    Test Sweep
     *************************/
/*
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
*/
    /********************************
     //    Residual
     *************************/
    double scale = 1/(h*h);
    DTMutableDoubleArray r(b.m(),b.n());
    double res = residual(u, b, scale, r);

    /********************************
     //    Test Residual
     *************************/
    DTMatlabDataFile outputResidual("OutputResidual.mat", DTFile::NewReadWrite);
      Write(outputResidual, "b", b);
      Write(outputResidual, "r", r);
      Write(outputResidual, "u", u);


    /********************************
     //    Coarsen
     *************************/
    int dimX = b.m()/2 + 1;
    int dimY = b.n()/2 + 1;
    DTMutableMesh2D coarseGrid(DTMesh2DGrid(origin, 2*h, 2*h, dimX, dimY), DTMutableDoubleArray(dimX,dimY));
    DTMesh2D fineGrid(DTMesh2DGrid(origin, h, h, b.m(), b.n()), r);
    coarsen(fineGrid, coarseGrid);

    /********************************
     //    Solve
     *************************/
    DTMutableDoubleArray rc = coarseGrid.DoubleData();
    DTMutableDoubleArray uc = DTMutableDoubleArray(rc.m(), rc.n());
    double* ucPointer = uc.Pointer();
    for(int i = 0; i < uc.Length(); ++i)
    {
        ucPointer[i] = 0;
    }
    recurse(uc, rc, omega, 2*h, Nb, Na, origin);

    /********************************
     //    Refine & Update
     *************************/
    DTMutableMesh2D targetCoarseGrid(DTMesh2DGrid(origin, 2*h, 2*h, uc.m(), uc.n()), uc);
    DTMutableMesh2D refineGrid(DTMesh2DGrid(origin, h , h, b.m(), b.n()), DTMutableDoubleArray(b.m(),b.n()));
    refine(targetCoarseGrid, refineGrid);
    DTMutableDoubleArray refine_residual = refineGrid.DoubleData();
    double* uPointer = u.Pointer();
    double* rPointer = refine_residual.Pointer();
    for(int i = 1; i < refineGrid.m() - 1; ++i)
    {
        for(int j = 1; j < refineGrid.n() - 1; ++j)
        {
            int idx = i*refineGrid.n() +j;
            u(i,j) += refine_residual(i,j);
//            uPointer[idx] += rPointer[idx];
        }
    }
    res = residual(u, b, scale, r);
    /********************************
     //    Sweep after
     *************************/
    sweep(u, b, omega, h, Na);
    res = residual(u, b, scale, r);

}
int main(int argc, char** argv)
{
    /********************************
     //    Read input data
     *************************/
    int Nb, Na, NList, id;
    double omega;
    if(argc < 6)
    {
        std::cout << "Usage: " << argv[0] << " <Length of NList>  <sweep before> <sweeps after> <omega> <instance id>" << std::endl;
        std::cout << "Currently use default value." << std::endl;
        Nb =  3;
        Na = 3;
        omega = 0.8;
        id = 1111;
    }
    else
    {
        NList = atoi(argv[1]);
        Nb = atoi(argv[2]);
        Na = atoi(argv[3]);
        omega = atof(argv[4]);
        id = atoi(argv[5]);

    }


    // Successive over relaxation(w) method
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
//    directSolve(b.m()-2, b.n()-2, h, b, u);

    double inf_norm;
//    inf_norm = residual(u, b.Pointer(), scale, r.Pointer());

    for(int i = 0; i < 100; ++i)
    {
        recurse(u, b, omega, h, Nb, Na, grid.Origin());
        inf_norm = residual(u, b, scale, r);
        residuals.push_back(inf_norm);

    }

    DTMutableDoubleArray rList(residuals.size(), 1);
    ostringstream plotId;
    plotId << "r";
    plotId << id;

    // ostringstream uId;
    // uId << "u";
    // uId << id;
    std::string sId = plotId.str();
//    std::string suId = uId.str();
    for(int i = 0; i < residuals.size(); ++i)
    {
        rList(i) = residuals[i];
    }

    Write(outputFile, sId, rList);
    Write(outputFile, "u", u);

    return 1;
}
