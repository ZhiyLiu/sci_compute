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

typedef Eigen::SparseMatrix<double> SpMat; 
typedef Eigen::Triplet<double> T;
void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);

int main(int argc, char** argv)
{
    //DTSetArguments(argc,argv);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    DTMatlabDataFile inputFile("Input.mat", DTFile::ReadOnly);
    DTMesh2D f;
    Read(inputFile, "f", f);
    DTFunction2D g;
    Read(inputFile, "g", g);

    DTMesh2DGrid grid = f.Grid();
    double h = grid.dx();
    DTDoubleArray fData = f.DoubleData();
    int cols = grid.n()-2;
    int n = (grid.m()-2)*(grid.n()-2);  // size of the image
    int gridWidth = grid.n();
    int gridHeight = grid.m();
    //test
    n = 25; gridWidth = 5; gridHeight = 5; DTMutableDoubleArray ab(n, n); DTMutableDoubleArray bt(n,1);
    // Assembly:
    std::vector<T> coefficients;            // list of non-zeros coefficients
    Eigen::VectorXd b(n);                   // the right hand side-vector resulting from the constraints
    double xzero = grid.Origin().x;
    double yzero = grid.Origin().y;

    // iterate over the grid, each element in the grid relates one row in coeff matrix
    for(int i = 0; i < gridHeight; ++i)
    {
        for(int j = 0; j < gridWidth; ++j)
        {
            // row number in coeff matrix
            int rowNum = i*gridWidth+j;
            double ytemp = yzero+j*h;
            double xtemp = xzero + i*h;
            bt(rowNum) = 2.4;
            b(rowNum) = 2.4;
            //for boundary coeff is 1
            if(i == 0 || i == gridHeight-1 || j==0 || j==gridWidth-1)
            {
                coefficients.push_back(T(rowNum, rowNum, 1));
                ab(rowNum, rowNum) = 1;
                b(rowNum) = rowNum;
                bt(rowNum) = rowNum;
            }
            else
            {
                coefficients.push_back(T(rowNum, rowNum, 4));
                ab(rowNum, rowNum) = 4;
                // right to diag
                coefficients.push_back(T(rowNum, rowNum+1,-1));
                coefficients.push_back(T(rowNum, rowNum+gridWidth, -1));
//                coefficients.push_back(T(rowNum+1, rowNum,-1));
//                coefficients.push_back(T(rowNum+gridWidth, rowNum, -1));
                ab(rowNum, rowNum+1) = -1;
                ab(rowNum, rowNum + gridWidth) = -1;

                // left to diag
                coefficients.push_back(T(rowNum, rowNum-1, -1));
                coefficients.push_back(T(rowNum, rowNum-gridWidth, -1));
//                coefficients.push_back(T(rowNum-1, rowNum, -1));
//                coefficients.push_back(T(rowNum-gridWidth, rowNum, -1));
                ab(rowNum, rowNum-1) = -1;
                ab(rowNum, rowNum-gridWidth) = -1;

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
    for(int i = 0; i < n; ++i)
    {
        std::cout << b(i) << " ";
    }
    std::cout << "x===" << endl;
    for(int i = 0; i < n; ++i)
    {
        std::cout << x(i) << " ";
    }

    std::cout << std::endl;
    std::cout << A << std::endl;
/*    for(int i = 0; i < n; ++i)
    {
        for(int j =0; j < n; ++j)
        {
            std::cout << A.coeff(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << endl;
*/
    //output
    DTMatlabDataFile outputFile("Output.mat", DTFile::NewReadWrite);
//    DTMesh2D uMesh(grid, returnArray);
//    Write(outputFile, "u", uMesh);
    outputFile.Save(ab, "Atest");
    outputFile.Save(bt, "btest");
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