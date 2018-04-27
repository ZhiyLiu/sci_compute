#include <math.h>
#include <iostream>
#include <algorithm>
#include "DTDoubleArray.h"
#include "DTMesh2D.h"

void coarsen(DTMesh2D& input, DTMutableMesh2D& output)
{
    DTDoubleArray inputData = input.DoubleData();//fine grid

    DTMutableDoubleArray outputData = output.DoubleData(); // coarse grid

    const double* inputPointer = inputData.Pointer();
    double* outputPointer = outputData.Pointer();
    for(int j= 0; j < inputData.n(); ++j)
    {
        for(int  i = 0; i < inputData.m(); ++i)
        {
            // keep the even row and column
            if(i % 2 == 0 && j % 2 == 0)
            {
                int inputIndex = i * inputData.n() + j;
                int iWidth = inputData.n();
                int outputIndex = 0.5 * i * outputData.n() + 0.5 * j;
                // boundary
                if(i == 0 || i==inputData.m()-1 || j == 0 || j == inputData.n() - 1)
                {
                    //outputData(i/2, j/2) = inputData(i,j);
                    outputPointer[outputIndex] = inputPointer[inputIndex];
                }
                else
                {
                    /* double sum = inputData(i,j) +\ */
                    /*     0.5 * inputData(i-1,j) +0.5 * inputData(i+1,j) + \ */
                    /*     0.5 * inputData(i,j-1) + 0.5 * inputData(i,j+1) +   \ */
                    /*     0.25 * inputData(i-1,j-1) + 0.25 * inputData(i+1,j+1) + 0.25 * inputData(i+1,j-1) +0.25 * inputData(i-1,j+1); */
                    /* outputData(i/2,j/2) = 0.25 * sum; */
                    double sum = inputPointer[inputIndex] + \
                        0.5 * inputPointer[inputIndex-1] + 0.5 * inputPointer[inputIndex+1] + \
                        0.5 * inputPointer[inputIndex-iWidth] + 0.5 * inputPointer[inputIndex+iWidth] + \
                        0.25 * inputPointer[inputIndex-iWidth-1] + 0.25 * inputPointer[inputIndex-iWidth+1] +\
                        0.25 * inputPointer[inputIndex+iWidth-1] + 0.25 * inputPointer[inputIndex+iWidth+1];
                    outputPointer[outputIndex] = 0.25 * sum;
                }

            }
        }
    }
}

void refine(DTMesh2D& input, DTMutableMesh2D& output)
{
    DTDoubleArray inputData = input.DoubleData();

    DTMutableDoubleArray outputData = output.DoubleData();
    const double* inputPointer = inputData.Pointer();
    double* outputPointer = outputData.Pointer();
    // reset value
    for(int j = 0; j < outputData.n(); ++j)
    {
        for(int  i = 0; i < outputData.m(); ++i)
        {
            int outputIndex = i * outputData.n() + j;
//            outputData(i,j) = 0.0;
            outputPointer[outputIndex] = 0.0;
        }
    }

    for(int  j = 0; j < outputData.n(); ++j)
    {
        for(int i = 0; i < outputData.m(); ++i)
        {
            if(i % 2 == 0 && j % 2 == 0)
            {
                int inputIndex = 0.5 * i * inputData.n() + 0.5 * j;
                int outputIndex = i * outputData.n() + j;
                int oWidth = outputData.n();

                double val = inputPointer[inputIndex];
                outputPointer[outputIndex] = val;
//                outputData(i,j) = inputData(i/2,j/2);
            }
        }
    }
    for(int  j = 0; j < outputData.n()-1; ++j)
    {
        for(int i = 0; i < outputData.m(); ++i)
        {
            if(i % 2 == 0 && j % 2 == 0)
            {
                int inputIndex = 0.5 * i * inputData.n() + 0.5 * j;
                int outputIndex = i * outputData.n() + j;
                int oWidth = outputData.n();
                outputPointer[outputIndex+1] = 0.5 * (inputPointer[inputIndex] + inputPointer[inputIndex+1]);
//                outputData(i,j+1) = 0.5 * (inputData(i/2,j/2) + inputData(i/2, j/2+1));
            }
        }
    }
    for(int  j = 0; j < outputData.n(); ++j)
    {
        for(int i = 0; i < outputData.m() - 2; ++i)
        {
            if(i % 2 == 0 && j % 2 == 0)
            {
                int inputIndex = 0.5 * i * inputData.n() + 0.5 * j;
                int outputIndex = i * outputData.n() + j;
                int oWidth = outputData.n();
                outputPointer[outputIndex + oWidth] = (outputPointer[outputIndex] + outputPointer[outputIndex+2*oWidth]) * 0.5;
//                outputData(i+1,j) = 0.5 * (inputData(i/2,j/2) + inputData(i/2+1, j/2));

            }
        }
    }
    for(int  j = 0; j < outputData.n() - 1; ++j)
    {
        for(int i = 0; i < outputData.m() - 1; ++i)
        {
            if(i % 2 == 0 && j % 2 == 0)
            {
                int inputIndex = 0.5 * i * inputData.n() + 0.5 * j;
                int outputIndex = i * outputData.n() + j;
                int oWidth = outputData.n();
                outputPointer[outputIndex + oWidth + 1] = 0.25 * (inputPointer[inputIndex] + inputPointer[inputIndex + 1] + inputPointer[inputIndex+inputData.n()] + inputPointer[inputIndex+inputData.n() + 1]);
//                outputData(i+1, j+1) = 0.25 * (inputData(i/2,j/2) + inputData(i/2+1,j/2) + inputData(i/2, j/2+1) + inputData(i/2+1, j/2+1));
            }
        }
    }
}
