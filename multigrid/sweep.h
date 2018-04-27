#include "DTDoubleArray.h"

// Weighted Jacobi
int sweep(DTMutableDoubleArray& u, DTMutableDoubleArray& b, double omega, double h, int N)
{
    double t = 0;

    double hsq = h * h;
    double qomega = 0.25 * omega;
    double oneMinusOmega = 1-omega;

    double* rhs = b.Pointer();
    DTMutableDoubleArray uNew = u.Copy();

    while(t < N)
    {
        t++;
        double* data = u.Pointer();
        double* dataNew = uNew.Pointer();
        for(int j = 1; j < u.n()-1; ++j)
        {
            for(int i = 1; i < u.m() -1; ++i)
            {
//                int idx = i * u.n() + j;
                int idx = j * u.n() + i;

//                uNew(i,j) = (1-omega) * u(i,j) + qomega * (u(i-1,j) + u(i+1,j) + u(i, j-1) + u(i,j+1) - hsq*b(i,j));

                /* double leftV = data[idx-1]; */
                /* double rightV = data[idx+1]; */
                /* double up = data[idx+u.n()]; */
                /* double d  = data[idx-u.n()]; */
                double v = (1-omega) * data[idx] + qomega * (data[idx -1] + data[idx+1] + data[idx + u.n()] + data[idx - u.n()] - hsq*rhs[idx]);
                dataNew[idx] = v;
            }
        }
        Swap(u, uNew);
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