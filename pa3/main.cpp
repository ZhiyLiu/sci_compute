#include <math.h>
#include <iostream>
#include "lapacke.h"
#include "DTMatlabDataFile.h"
#include "DTArguments.h"
#include "DTDoubleArray.h"
#include "DTIntArray.h"

int main(int argc, const char *argv[])
{
    DTSetArguments(argc,argv);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    DTMatlabDataFile inputFile("Input.mat", DTFile::ReadOnly);
    DTDoubleArray xy = inputFile.ReadDoubleArray("xy");
    DTDoubleArray x = inputFile.ReadDoubleArray("x");

    if(xy.m() != 2 || x.m() != 1 || xy.n() != 5)
    {
        std::cout << "Invalid input array dimensions." << std::endl;
        return -1;
    }

    // suppose the sample interval is the same
    double h = fabs(xy(0, 1) - xy(0,0));
    double hsq = h * h;
    int N = xy.n();
    double DL[N-1], DU[N-1], D[N], B[N];
    int INFO = 0;
    int NRHS = 1;
    int LDB = N;

    for(int i = 0; i < N-2; ++i)
    {
        DL[i] = h;
    }
    // the last elemt of DL is 0
    DL[N-2] = 0;

    // the first elemt of DU is 0
    DU[0] = 0;
    for(int j = 1; j < N-1; ++j)
    {
        DU[j] = h;
    }

    // the first and last are 1
    D[0] = D[N-1] = 1;
    for(int k = 1; k < N-1; ++k)
    {
        D[k] = 4 * h;
    }

    // the first  of B are 0
    B[0] = B[N-1] = 0;
    double ratio = 3/hsq;
    for(int m = 1; m < N-1; ++m)
    {
        double temp = xy(1,m);
        B[m] = ratio * (xy(1, m+1) - 2*xy(1, m) + xy(1, m-1));
    }

    LAPACK_dgtsv(&N, &NRHS, DL, D, DU, B, &LDB, &INFO);

    double coeffD[N-1], coeffB[N-1];
    ratio = 1/(3*h);
    for(int i = 0; i < N-1; ++i)
    {
        coeffD[i] = (B[i+1] - B[i]) * ratio;
        coeffB[i] = (xy(1, i + 1) - xy(1, i))/h - h*(B[i+1] + 2* B[i])/3;
    }

//    coeffB[N-1] = coeffB[N-2] + h * (B[N-2] + B[N-1]);
//    double hcub = h * h * h;
    //  coeffD[N-1] = (xy(1, N-1) - xy(1, N-2) - coeffB[N-1] * h) / hcub;

    int cols = x.n();
    DTMutableDoubleArray outputMat(1, cols);
    double *y = outputMat.Pointer();
    for(int i = 0; i < x.n(); ++i)
    {
        double temp = x(i);
        if(x(i) < xy(0,0))
        {
            // out of bound
            double dist = x(i) - xy(0,0);
            // s0
            y[i] = xy(1, 0) + coeffB[0] * dist + B[0] * dist * dist + coeffD[0] * dist * dist * dist;
        }
        if(x(i) >= xy(0, 0) && x(i) < xy(0, 1))
        {
            double dist = x(i) - xy(0,0);
            //s0
            y[i] = xy(1, 0) + coeffB[0] * dist + B[0] * dist * dist + coeffD[0] * dist * dist * dist;
        }

        if(x(i) >= xy(0, 1) && x(i) < xy(0, 2))
        {
            double dist = x(i) - xy(0,1);
            //s1
            y[i] = xy(1, 1) + coeffB[1] * dist + B[1] * dist * dist + coeffD[1] * dist * dist * dist;
        }

        if(x(i) >= xy(0, 2) && x(i) < xy(0,3))
        {
            double dist = x(i) - xy(0,2);
            //s2
            y[i] = xy(1, 2) + coeffB[2] * dist + B[2] * dist * dist + coeffD[2] * dist * dist * dist;
        }

        if(x(i) >= xy(0, 3) && x(i) < xy(0,4))
        {
            double dist = x(i) - xy(0,3);
            //s3
            y[i] = xy(1, 3) + coeffB[3] * dist + B[3] * dist * dist + coeffD[3] * dist * dist * dist;
        }

        if(x(i) >= xy(0,4))
        {
            double dist = x(i) - xy(0,4);
            // s3
            y[i] = xy(1, 4) + coeffB[3] * dist + B[3] * dist * dist + coeffD[3] * dist * dist * dist;
        }
    }

    DTMatlabDataFile outFile("Output.mat",DTFile::NewReadWrite);
    outFile.Save(outputMat, "y");
    std::cout << "program terminated."  << std::endl;

    return 0;
}