#include "DTDoubleArray.h"

// Weighted Jacobi
int sweep(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double omega, double h, int N)
{
    double t = 0;

    double hsq = h * h;
    double qomega = 0.25 * omega;
    double oneMinusOmega = 1-omega;
    double* data = u.Pointer();
    double* rhs = b.Pointer();
    while(t < N)
    {
        t++;
        for(int i = 1; i < u.m()-1; ++i)
        {
            for(int j = 1; j < u.n()-1; ++j)
            {
                int idx = i*u.n() + j;
/*                int idxL = idx-1;
                int idxR = idx+1;
                int idxU = idx-u.n();
                int idxD = idx+u.n();
                double u0 = u(idx);
                double uL = u(idxL);
                double uR = u(idxR);
                double uU = u(idxU);
                double uD = u(idxD);
*/
                u(i,j) = (1-omega) * u(i,j) + qomega * (u(i-1,j) + u(i+1,j) + u(i, j-1) + u(i,j+1) - hsq*b(i,j));

//                u(idx) = oneMinusOmega * u(idx) + qomega * (u(idx+1) + u(idx-1) + u(idx+u.n()) + u(idx-u.n()) - hsq * b(idx));
//                data[idx] = oneMinusOmega * data[idx] + qomega * (data[idx+1] + data[idx-1] + data[idx+u.n()] + data[idx-u.n()] - hsq * rhs[idx]);
            }
        }
    }
    return 1;
}

// Red-black Gauss Seidel
int sweep(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double omega, double h, int N, DTMesh2DGrid& grid)
{
    // read in reference solution from sparse solver

    DTMatlabDataFile referFile("ExactSolution.mat", DTFile::ReadOnly);
    DTMesh2D xexact;
    Read(referFile, "u", xexact);
    DTMutableDoubleArray x_exact = xexact.DoubleData().Copy();

    DTMatlabDataFile outputFile("Output.mat",DTFile::NewReadWrite);
    DTSeriesMesh2D computed(outputFile,"MyVar");

    double t = 0;
    computed.Add(DTMesh2D(grid,u),t); // Saves the time value to disk

    // Iterate, and increment t
    DTMatlabDataFile outputErrorFile("Errors.mat",DTFile::NewReadWrite);
    ////////error file intitial error 
    DTMutableDoubleArray matError(N+1, 1);
    double initError = -1000.0;
    for(int i = 0; i < u.m(); ++i)
    {
        for(int j = 0; j < u.n(); ++j)
        {
            double err = abs(u(i,j) - x_exact(i,j));
            if(err > initError)
            {
                initError = err;
            }
        }
    }

    matError(0) = initError;

    double hsq = h*h;

    while(t < N)
    {
        t += 1;
        double max_error = -1000.0;
        for(int i = 1; i < u.m()-1; ++i)
        {
            for(int j = 1; j < u.n()-1; ++j)
            {
                // sum is odd
                if((i+j)%2 != 0)
                {
                    double xnew = (1-omega)*u(i,j) +0.25* omega*(u(i-1,j) + u(i+1,j)+u(i,j-1)+u(i,j+1) - hsq*b(i,j));
                    u(i,j) = xnew;
                    double err = abs(xnew - x_exact(i,j));
                    if(err > max_error)
                    {
                        max_error = err;
                    }
                }
            }
        }
        for(int i = 1; i < u.m()-1; ++i)
        {
            for(int j = 1; j < u.n()-1; ++j)
            {
                // sum is even
                if((i+j)%2==0)
                {
                    double xnew = (1-omega)*u(i,j) + 0.25* omega*(u(i-1,j) + u(i+1,j)+u(i,j-1)+u(i,j+1) - hsq*b(i,j));
                    u(i,j) = xnew;
                    double err = abs(xnew - x_exact(i,j));
                    if(err > max_error)
                    {
                        max_error = err;
                    }
                }
            }
        }

        computed.Add(DTMesh2D(grid,u),t);
        matError((int)t) = max_error;
    }


    outputErrorFile.Save(matError, "errors");
    outputFile.Save(matError, "ExecutionErrors");
    return 0;
}