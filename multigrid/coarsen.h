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
#include "DTSeriesMesh2D.h"
#include "time.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

void coarsen(DTMesh2D& input, DTMutableMesh2D& output)
{
    DTDoubleArray inputData = input.DoubleData();//fine grid

    DTMutableDoubleArray outputData = output.DoubleData(); // coarse grid

    for(int i = 0; i < inputData.m(); ++i)
    {
        for(int j= 0; j < inputData.n(); ++j)
        {
            // keep the even row and column
            if(i % 2 == 0 && j % 2 == 0)
            {
                // boundary
                if(i == 0 || i==inputData.m()-1 || j == 0 || j == inputData.n() - 1)
                {
                    outputData(i/2, j/2) = inputData(i,j);
                }
                else
                {
                    double sum = inputData(i,j) +\
                        0.5 * inputData(i-1,j) +0.5 * inputData(i+1,j) + \
                        0.5 * inputData(i,j-1) + 0.5 * inputData(i,j+1) +   \
                        0.25 * inputData(i-1,j-1) + 0.25 * inputData(i+1,j+1) + 0.25 * inputData(i+1,j-1) +0.25 * inputData(i-1,j+1);
                    outputData(i/2,j/2) = 0.25 * sum;
                }

            }
        }
    }
}

void refine(DTMesh2D& input, DTMutableMesh2D& output)
{
    DTDoubleArray inputData = input.DoubleData();

    DTMutableDoubleArray outputData = output.DoubleData();

    // reset value
    for(int i = 0; i < outputData.m(); ++i)
    {
        for(int j = 0; j < outputData.n(); ++j)
        {
            outputData(i,j) = 0.0;
        }
    }
    for(int i = 1; i < outputData.m()-1; ++i)
    {
        for(int j= 1; j < outputData.n()-1; ++j)
        {
            if(i % 2 == 0 && j % 2 == 0)
            {
                outputData(i,j) = inputData(i/2, j/2);
                // distribute the weighted value to adjacent grid points
                outputData(i-1, j) += 0.5 * inputData(i/2, j/2);
                outputData(i+1, j) += 0.5 * inputData(i/2, j/2);
                outputData(i, j-1) += 0.5 * inputData(i/2, j/2);
                outputData(i, j+1) += 0.5 * inputData(i/2, j/2);

                outputData(i-1, j - 1) += 0.25 * inputData(i/2, j/2);
                outputData(i+1, j - 1) += 0.25 * inputData(i/2, j/2);
                outputData(i-1, j + 1) += 0.25 * inputData(i/2, j/2);
                outputData(i+1, j + 1) += 0.25 * inputData(i/2, j/2);
            }

        }
    }

}
