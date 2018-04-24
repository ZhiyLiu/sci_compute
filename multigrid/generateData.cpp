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
#include "DTPoint2D.h"
#include "time.h"

int main(int argc, char** argv)
{

    int dim = atoi(argv[1]);

    double h = 0.025;
    DTMesh2DGrid grid(DTPoint2D(0.0,0.0),h,h,dim,dim);

    DTMutableDoubleArray data(dim, dim);
    for(int i = 0; i < dim; ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            double x = i* h;
            double y = j * h;
            data(i,j) = sin((x *x + y*y)/dim);
        }
    }

    DTMesh2D f(grid, data);
    DTMatlabDataFile outputFile("MyInput.mat", DTFile::NewReadWrite);
    Write(outputFile, "f", f);

    return 1;
}