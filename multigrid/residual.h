#include <math.h>
#include <iostream>
#include <algorithm>
#include "DTDoubleArray.h"

//void construct_A(int m, int n, double scale, int cols);
double residual(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double scale, DTMutableDoubleArray& r)
{
    double inf_norm = -1.0;
    double* data = u.Pointer();
    for(int i = 1; i < u.m()-1; ++i)
    {
        for(int j = 1; j < u.n()-1; ++j)
        {
            // rowNum in right hand side vector
            int idx = i*u.n() + j;
            double ax = 4 * u(i,j) - u(i,j-1) - u(i,j+1) - u(i+1,j) - u(i-1,j);
//            double ax = 4 * data[idx] - data[idx -1] - data[idx + 1] - data[idx+u.n()] - data[idx - u.n()];
            ax *= scale;
            r(i,j) = -b(i,j) - ax;
            if(inf_norm < abs(r(i,j)))
            {
                inf_norm = abs(r(i,j));
            }
        }
    }
    return inf_norm;
}
// return inf norm
/* double residual(DTMutableDoubleArray& u, double* b, double scale, double* r) */
/* { */
/*     double inf_norm = -1.0; */
/*     double* data = u.Pointer(); */
/*     for(int i = 1; i < u.m()-1; ++i) */
/*     { */
/*         for(int j = 1; j < u.n()-1; ++j) */
/*         { */
/*             // rowNum in right hand side vector */
/*             int idx = i*u.n() + j; */
/*             double ax = 4 * u(i,j) - u(i,j-1) - u(i,j+1) - u(i+1,j) - u(i-1,j); */
/* //            double ax = 4 * data[idx] - data[idx -1] - data[idx + 1] - data[idx+u.n()] - data[idx - u.n()]; */
/*             ax *= scale; */
/*             r[idx] = -b[idx] - ax; */
/*             if(inf_norm < abs(r[idx])) */
/*             { */
/*                 inf_norm = abs(r[idx]); */
/*             } */
/*         } */
/*     } */
/*     return inf_norm; */
/* /\* */
/*     if(needTest) */
/*     { */
/*         construct_A(r.m(), r.m(), scale, u.n()); */
/*         DTMatlabDataFile outputFile("residual.mat", DTFile::NewReadWrite); */
/*         outputFile.Save(r, "r"); */
/*     } */
/* *\/ */
/* } */

void construct_A(int m, int n, double scale, int cols)
{
    DTMutableDoubleArray a(m,n);
    for(int j = 0; j < n;++j)
    {
        for (int i = 0; i < m; ++i)
        {
            if(i == j)
            {
                a(i,j) = 4 * scale;
            }
            if(j-i == 1)
            {
                if((i+1)%(cols) == 0)
                {
                    a(i,j) = 0;
                }
                else{
                    a(i,j) = -scale;
                }

            }
            if(j-i > 1 && j-i < cols)
            {
                a(i,j) = 0;
            }
            if(j-i == cols)
            {
                a(i,j) = -scale;
            }
            if(i-j == 1)
            {
                if((j+1)%(cols)==0)
                {
                    a(i,j) = 0;
                }
                else{
                    a(i,j) = -scale;
                }
            }
            if(i-j >1 && i-j < cols)
            {
                a(i,j) = 0;
            }
            if(i-j == cols)
            {
                a(i,j) = -scale;
            }
        }
    }
    DTMatlabDataFile outputFile("Operater.mat",DTFile::NewReadWrite);
    outputFile.Save(a, "A");
}