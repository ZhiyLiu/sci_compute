#include <math.h>
#include <iostream>
#include <algorithm>
#include "DTDoubleArray.h"
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);

// gridHeight & gridWidth don't include bdry
// fData & u include bdry
void directSolve(int gridHeight, int gridWidth, double h, DTDoubleArray& fData, DTMutableDoubleArray& u)
{
    std::vector<T> coefficients;            // list of non-zeros coefficients
    int n = gridWidth * gridHeight;
    Eigen::VectorXd b(n);                   // the right hand side-vector resulting from the constraints
    double hsq = h* h;
    for(int i = 0; i < gridHeight; ++i)
    {
        for(int j = 0; j < gridWidth; ++j)
        {
            // row number in coeff matrix
            int rowNum = i*gridWidth+j;

            coefficients.push_back(T(rowNum, rowNum, 4));
            //ab(rowNum, rowNum) = 4;
            double temp = fData(i+1, j+1);
            b(rowNum) = 0 - fData(i+1,j+1) * hsq;
            //          bt(rowNum) = b(rowNum);

            // one side
            if(i==0)
            {
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));

                if(j ==0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //    ab(rowNum, rowNum+1) = -1;
                }
                else if(j == gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;
                }
                else
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    // ab(rowNum, rowNum-1) = -1;
                    // ab(rowNum, rowNum+1) = -1;
                }

            }
            else if(i == gridHeight-1)
            {
                coefficients.push_back(T(rowNum, rowNum - gridWidth, -1));
                //ab(rowNum, rowNum - gridWidth) = -1;

                if(j == 0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //  ab(rowNum, rowNum+1) = -1;

                }
                else if(j==gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;
                }
                else
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    coefficients.push_back(T(rowNum, rowNum+1, -1));

                    // ab(rowNum, rowNum-1) = -1;
                    // ab(rowNum, rowNum+1) = -1;
                }
            }
            else if(j == 0)
            {
                // corners have been tackled above
                coefficients.push_back(T(rowNum, rowNum+1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum+1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;

            }
            else if(j == gridWidth-1)
            {
                coefficients.push_back(T(rowNum, rowNum -1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum-1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;

            }
            else
            {
                coefficients.push_back(T(rowNum, rowNum+1, -1));
                coefficients.push_back(T(rowNum, rowNum -1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum+1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;
                // ab(rowNum, rowNum-1) = -1;
            }
        }
    }

    SpMat A(n,n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solving:
    Eigen::SimplicialLLT<SpMat> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success)
    {
        return;
    }
    Eigen::VectorXd x = solver.solve(b);         // use the factorization to solve for the given right hand side

    for(int i = 0; i < x.size(); ++i)
    {
        int r = i / gridWidth;
        int c = i % gridWidth;
        u(r+1,c+1) = x[i];
    }

}
int sparseSolve()
{
    //DTSetArguments(argc,argv);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    DTMatlabDataFile inputFile("MyInput.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);

    DTMesh2DGrid grid = f.Grid();
    double h = grid.dx();
    DTDoubleArray fData = f.DoubleData();
    int cols = grid.n()-2;
    int n = (grid.m()-2)*(grid.n()-2);  // size of the image
    int gridWidth = grid.n()-2;
    int gridHeight = grid.m()-2;
    //test
    //n = 16; gridWidth = 4; gridHeight = 4;
//    DTMutableDoubleArray ab(n, n); DTMutableDoubleArray bt(n,1);
    DTMutableDoubleArray returnArray(gridHeight+2, gridWidth+2);
    // Assembly:
    std::vector<T> coefficients;            // list of non-zeros coefficients
    Eigen::VectorXd b(n);                   // the right hand side-vector resulting from the constraints
    double xzero = grid.Origin().x;
    double yzero = grid.Origin().y;
    double hsq = h * h;
    // iterate over the grid, each element in the grid relates one row in coeff matrix
    for(int i = 0; i < gridHeight; ++i)
    {
        for(int j = 0; j < gridWidth; ++j)
        {
            // row number in coeff matrix
            int rowNum = i*gridWidth+j;

            coefficients.push_back(T(rowNum, rowNum, 4));
            //ab(rowNum, rowNum) = 4;
            double temp = fData(i+1, j+1);
            b(rowNum) = 0 - fData(i+1,j+1) * hsq;
            //          bt(rowNum) = b(rowNum);

            // one side
            if(i==0)
            {
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));

                if(j ==0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //    ab(rowNum, rowNum+1) = -1;
                }
                else if(j == gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;
                }
                else
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    // ab(rowNum, rowNum-1) = -1;
                    // ab(rowNum, rowNum+1) = -1;
                }

            }
            else if(i == gridHeight-1)
            {
                coefficients.push_back(T(rowNum, rowNum - gridWidth, -1));
                //ab(rowNum, rowNum - gridWidth) = -1;

                if(j == 0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //  ab(rowNum, rowNum+1) = -1;

                }
                else if(j==gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;
                }
                else
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    coefficients.push_back(T(rowNum, rowNum+1, -1));

                    // ab(rowNum, rowNum-1) = -1;
                    // ab(rowNum, rowNum+1) = -1;
                }
            }
            else if(j == 0)
            {
                // corners have been tackled above
                coefficients.push_back(T(rowNum, rowNum+1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum+1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;

            }
            else if(j == gridWidth-1)
            {
                coefficients.push_back(T(rowNum, rowNum -1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum-1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;

            }
            else
            {
                coefficients.push_back(T(rowNum, rowNum+1, -1));
                coefficients.push_back(T(rowNum, rowNum -1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum+1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;
                // ab(rowNum, rowNum-1) = -1;
            }
        }
    }

    SpMat A(n,n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solving:
    Eigen::SimplicialLLT<SpMat> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success)
    {
        return -1;
    }
    Eigen::VectorXd x = solver.solve(b);         // use the factorization to solve for the given right hand side


    // corners
    for(int i = 0; i < x.size(); ++i)
    {
        int r = i / gridWidth;
        int c = i % gridWidth;
        returnArray(r+1,c+1) = x[i];
    }

    //output
    DTMatlabDataFile outputFile("ExactSolution.mat", DTFile::NewReadWrite);
    DTMesh2D uMesh(grid, returnArray);
    Write(outputFile, "u", uMesh);
//    outputFile.Save(ab, "Atest");
//    outputFile.Save(bt, "btest");
    return 0;
}
void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs,
                       Eigen::VectorXd& b, const Eigen::VectorXd& boundary)
{
    int n = int(boundary.size());
    int id1 = i+j*n;
    if(i==-1 || i==n) b(id) -= w * boundary(j); // constrained coefficient
    else  if(j==-1 || j==n) b(id) -= w * boundary(i); // constrained coefficient
    else  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}
void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n)
{
    b.setZero();
    Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);
    for(int j=0; j<n; ++j)
    {
        for(int i=0; i<n; ++i)
        {
            int id = i+j*n;
            insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
            insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
            insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
            insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
            insertCoefficient(id, i,j,    4, coefficients, b, boundary);
        }
    }
}