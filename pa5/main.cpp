#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <vector>
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
#include "time.h"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);

int main(int argc, char** argv)
{
    //DTSetArguments(argc,argv);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    DTMatlabDataFile inputFile("MyInput.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);
    DTFunction2D g;
    Read(inputFile, "g", g);

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
                //ab(rowNum, rowNum+gridWidth) = -1;
                double xtemp = xzero;
                double ytemp = yzero+(j+1)*h;
                b(rowNum) += g(xtemp, ytemp);
                //  bt(rowNum) += g(xtemp, ytemp);
                returnArray(i,j+1) = g(xtemp, ytemp);
                if(j ==0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //    ab(rowNum, rowNum+1) = -1;
                    b(rowNum) += g(xzero+h, yzero);
                    //  bt(rowNum) += g(xzero+h, yzero);
                    returnArray(i+1, j) = g(xzero+h, yzero);
                }
                else if(j == gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;

                    b(rowNum) += g(xzero+h, ytemp);
                    //bt(rowNum) += g(xzero+h, ytemp);
                    returnArray(i+1, j+2) = g(xzero+h, ytemp);
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

                double xtemp = xzero + (gridHeight+1)*h;
                double ytemp = yzero + (j+1)*h;
                b(rowNum) += g(xtemp, ytemp);
                //bt(rowNum) +=g(xtemp, ytemp);
                returnArray(i+2,j+1) = g(xtemp,ytemp);
                if(j == 0)
                {
                    coefficients.push_back(T(rowNum, rowNum+1, -1));
                    //  ab(rowNum, rowNum+1) = -1;

                    b(rowNum) += g(xtemp-h, yzero);
                    //  bt(rowNum) += g(xtemp-h, yzero);

                    returnArray(i+1, 0) = g(xtemp-h, yzero);
                }
                else if(j==gridWidth-1)
                {
                    coefficients.push_back(T(rowNum, rowNum-1, -1));
                    //ab(rowNum, rowNum-1) = -1;

                    b(rowNum) += g(xtemp-h, ytemp+h);
                    //bt(rowNum) += g(xtemp-h, ytemp+h);
                    returnArray(i+1, j+2) = g(xtemp-h,ytemp+h);
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

                double xtemp = xzero + (i+1)*h;
                double ytemp = yzero;
                b(rowNum) += g(xtemp, ytemp);
                //bt(rowNum) += g(xtemp, ytemp);

                returnArray(i+1,j)=g(xtemp,ytemp);
            }
            else if(j == gridWidth-1)
            {
                coefficients.push_back(T(rowNum, rowNum -1, -1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
                // ab(rowNum, rowNum-gridWidth) = -1;
                // ab(rowNum, rowNum-1) = -1;
                // ab(rowNum, rowNum+gridWidth) = -1;

                double xtemp = xzero + (i+1)*h;
                double ytemp = yzero + (gridWidth+1)*h;
                b(rowNum) += g(xtemp, ytemp);
                //bt(rowNum) += g(xtemp, ytemp);

                returnArray(i+1, j+2) = g(xtemp,ytemp);
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
    clock_t t1,t2;
    t1 = clock();
    Eigen::SimplicialLLT<SpMat> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success)
    {
        return -1;
    }
    Eigen::VectorXd x = solver.solve(b);         // use the factorization to solve for the given right hand side
    t2 = clock();
    std::cout << "Time: " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    // corners
    int t = grid.n() -1;
    int q = grid.m() - 1;
    returnArray(0,0) = g(0,0);
    returnArray(0, t) = g(xzero,(grid.n()-1)*h+yzero);
    returnArray(q, 0) = g((grid.m()-1)*h+xzero, yzero);
    returnArray(q, t) = g((grid.m()-1)*h+xzero,(grid.n()-1)*h+yzero);

    for(int i = 0; i < x.size(); ++i)
    {
        int r = i / gridWidth;
        int c = i % gridWidth;
        returnArray(r+1,c+1) = x[i];
    }

    //output
    DTMatlabDataFile outputFile("Output.mat", DTFile::NewReadWrite);
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