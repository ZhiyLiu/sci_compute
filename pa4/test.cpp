#include <math.h>
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

int main(int argc, const char *argv[])
{

    DTMatlabDataFile inputFile("Input.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);
    DTFunction2D g;
    Read(inputFile, "g", g);

    DTMesh2DGrid grid = f.Grid();
    double h = grid.dx();
    DTDoubleArray fData = f.DoubleData();

//    DTFunction2D x = DTFunction2D::x();
//    DTFunction2D y = DTFunction2D::y();

//    DTFunction2D g = x*x+y;
    if(grid.n() < 3)
    {
        std::cout << "invalid parameters." << std::endl;
        return -1;
    }

    int kl = grid.n() - 2;
    int ku = grid.n() - 2;
    int nrhs = 1;
    int n = (grid.m() - 2) * (grid.n() - 2);
    int ldb = n;
    int ldab = 2*kl+ku +1;
    int rows = 2*kl + ku + 1;
    int cols = grid.n()-2;

    int ipiv[n];
    int info;

/*    int kl = 4, ku = 4, n = 16;
    int rows = 2*kl+ku+1;
    int cols = 4;
*/
    DTMutableDoubleArray AB(rows, n);
    for(int i = 0; i < n; i++)
    {
        AB(kl+ku, i) = 4;
        AB(kl, i) = -1;

        if(i % cols == 0 && i != n-1)
        {
            AB(kl+ku-1, i) = 0;
        }
        else if(i>0)
        {
            AB(kl+ku-1, i) = -1;
        }

        for(int j = kl+1; j < kl+ku-1; ++j)
        {
            AB(j, i) = 0;
        }

        if((i+1)%cols == 0 && i != 0)
        {
            AB(kl+ku+1,i) = 0;
        }
        else
        {
            AB(kl+ku+1, i) = -1;
        }

        for(int j = kl+ku+2; j < rows; ++j)
        {
            AB(j, i) = 0;
        }
        AB(rows-1,i) = -1;

    }
    DTMutableDoubleArray a(n,n);
    for(int j = 0; j < n;++j)
    {
        for (int i = 0; i < n; ++i)
        {
            if(i == j)
            {
                a(i,j) = 4;
            }
            if(j-i == 1)
            {
                if((i+1)%(cols) == 0)
                {
                    a(i,j) = 0;
                }
                else{
                    a(i,j) = -1;
                }

            }
            if(j-i > 1 && j-i < cols)
            {
                a(i,j) = 0;
            }
            if(j-i == cols)
            {
                a(i,j) = -1;
            }
            if(i-j == 1)
            {
                if((j+1)%(cols)==0)
                {
                    a(i,j) = 0;
                }
                else{
                    a(i,j) = -1;
                }
            }
            if(i-j >1 && i-j < cols)
            {
                a(i,j) = 0;
            }
            if(i-j == cols)
            {
                a(i,j) = -1;
            }
        }
    }

    DTMutableDoubleArray b(n);
    for(int i = 0; i < n; ++i)
    {
        b(i) = 2.4;
    }
    DTMatlabDataFile outputFile("Output.mat");
//    outputFile.Save(a, "An");
//    outputFile.Save(b,"bn");

    LAPACK_dgbsv(&n, &kl, &ku, &nrhs, AB.Pointer(), &ldab, ipiv, b.Pointer(), &ldb, &info);
    outputFile.Save(b, "xx");
    return 0;
}