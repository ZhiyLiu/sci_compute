#include "stdlib.h"
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "lapacke.h"
#include "DTMatlabDataFile.h"
#include "DTArguments.h"
#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTMesh2D.h"

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

std::vector<double> residuals;
void recurse(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double omega, double h, int Nb, int Na, DTPoint2D origin)
{
    // terminate criteria
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
    double inf_norm = residual(u, b.Pointer(), scale, r.Pointer());
    residuals.push_back(inf_norm);

    /********************************
     //    Test Residual
     *************************/
/*    DTMatlabDataFile outputResidual("OutputResidual.mat", DTFile::NewReadWrite);
      Write(outputResidual, "b", b);
      Write(outputResidual, "r", r);
      Write(outputResidual, "u", u);
*/

    /********************************
     //    Test Coarsen & Refine
     *************************/
/*    int dim = 9; // 9*9 include bdry should coarsen to 5*5
      h = 0.025;
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
*/

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
    recurse(uc, rc, omega, 2*h, Nb, Na, origin);

    /********************************
     //    Refine & Update
     *************************/
    DTMutableMesh2D refineGrid(DTMesh2DGrid(origin, h , h, b.m(), b.n()), DTMutableDoubleArray(b.m(),b.n()));
    refine(coarseGrid, refineGrid);
    DTMutableDoubleArray refine_residual = refineGrid.DoubleData();
    for(int i = 1; i < refineGrid.m() - 1; ++i)
    {
        for(int j = 1; j < refineGrid.n() - 1; ++j)
        {
            u(i,j) += refine_residual(i,j);
        }
    }

    /********************************
     //    Sweep after
     *************************/
    sweep(u, b, omega, h, Na);

}
int main(int argc, char** argv)
{
    /********************************
     //    Read input data
     *************************/
    if(argc < 6)
    {
        std::cout << "Usage: " << argv[0] << " <Length of NList>  <sweep before> <sweeps after> <omega> <#iter>" << std::endl;
        return -1;
    }

    int NList = atoi(argv[1]);
    int Nb = atoi(argv[2]);
    int Na = atoi(argv[3]);
    double omega = atof(argv[4]);
    int iterNum = atoi(argv[5]);
/*    int Nb =  3;
      int Na = 3;
      double omega = 0.8;
*/
    // Successive over relaxation(w) method
    DTMatlabDataFile inputFile("MyInput.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);
    DTMesh2D u0;
    Read(inputFile, "u", u0);

    DTMesh2DGrid grid = f.Grid();
    double h = grid.dx();
    DTMutableDoubleArray u = u0.DoubleData().Copy();
    DTMutableDoubleArray b = f.DoubleData().Copy();

    recurse(u, b, omega, h, Nb, Na, grid.Origin());

    DTMatlabDataFile outputFile("OutputAll.mat", DTFile::NewReadWrite);

    DTMutableDoubleArray rList(residuals.size());
    if(outputFile.Contains("r"))
    {
        Read(outputFile, "r", rList);
    }

    for(int i = 0; i < residuals.size(); ++i)
    {
        rList(i) = residuals[i];
    }
    Write(outputFile, "r", rList);
    Write(outputFile, "u", u);

    return 1;
}
