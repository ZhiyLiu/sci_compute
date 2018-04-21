#include <iostream>
#include <math.h>
#include <lapacke.h>
#include "DTMatlabDataFile.h"
#include "DTArguments.h"
#include "DTDoubleArray.h"
#include "DTIntArray.h"
char    TRANS = 'N';
int calculateFactorization(DTDoubleArray A, DTDoubleArray B)
{

    DTMutableDoubleArray mA = A.Copy();
    DTMutableDoubleArray mB = B.Copy();

    int M = mA.m();
    int N = mA.n();
    int NRHS = mB.n();
    int LDA = M;
    int LDB = mB.m();
    int minDim = min(M,N);
    DTMutableIntArray IPIV(minDim);
    int INFO = 3;
    LAPACK_dgetrf(&M,&N,mA.Pointer(),&LDA,IPIV.Pointer(),&INFO);


    DTMatlabDataFile outputFile("Output.mat",DTFile::NewReadWrite);

    // Output from LU & P
    outputFile.Save(mA,"LU");
    outputFile.Save(IPIV, "P");

    return 0;
}

int solveSystem(DTDoubleArray A, DTDoubleArray B)
{
    DTMutableDoubleArray mA = A.Copy();
    DTMutableDoubleArray mB = B.Copy();
    int INFO = 3;
    int M = mA.m();
    int N = mA.n();
    int NRHS = mB.n();
    int LDA = M;
    int LDB = mB.m();
    int minDim = min(M,N);
    DTMutableIntArray IPIV(minDim);
    
    LAPACK_dgetrf(&M,&N,mA.Pointer(),&LDA,IPIV.Pointer(),&INFO);
    LAPACK_dgetrs(&TRANS,&N,&NRHS,mA.Pointer(),&LDA,IPIV.Pointer(),mB.Pointer(),&LDB,&INFO);

    DTMatlabDataFile outFile2("Output.mat");
    outFile2.Save(mB, "x");

    return 0;
}

int estimateConditionNumber(DTDoubleArray A)
{

    // estimate condition number
    char NORM = '1'; // one-norm
    int INFO = 3;
    double ANORM = 0.0;
    // one norm of A
    for(int i = 0; i < A.n(); ++i)
    {
        double tempSum = 0.0;
        for(int j = 0; j < A.m(); ++j)
        {
            tempSum += fabs(A(i,j));
        }
        if(tempSum > ANORM)
        {
            ANORM = tempSum;
        }
    }
    double RCOND = 0.0;
    DTMutableDoubleArray mA = A.Copy();
    
    int M = mA.m();
    int N = mA.n();
    int LDA = M;
    int minDim = min(M,N);
    int IPIV[minDim];
    int IWORK[N];
    double WORK[4*N];

    LAPACK_dgetrf(&M,&N,mA.Pointer(),&LDA,IPIV,&INFO);

    LAPACK_dgecon(&NORM, &N, mA.Pointer(), &LDA, &ANORM, &RCOND, WORK, IWORK, &INFO );

    if(INFO < 0)
    {
        std::cout << "Illegal value at " << INFO << " th argument." << std::endl;
        return -1;
    }

    double INFNORM = 0.0;
    for(int i = 0; i < A.m(); ++i)
    {
        double tempSum = 0.0;
        for(int j = 0; j < A.n(); ++j)
        {
            tempSum += fabs(A(i,j));
        }
        if(tempSum > INFNORM)
        {
            INFNORM = tempSum;
        }
    }
    
    double RCOND2 = 0.0;
    int IIWORK[N];
    double oWork[4*N];
    char NORM2 = 'I';
    LAPACK_dgecon(&NORM2, &N, mA.Pointer(), &LDA, &INFNORM, &RCOND2, oWork, IIWORK, &INFO );
    if(INFO < 0)
    {
        std::cout << "Illegal value at " << INFO << " th argument." << std::endl;
        return -1;

    }

    DTMatlabDataFile outputFile("Output.mat");
    outputFile.Save(RCOND, "KappaOne");
    outputFile.Save(RCOND2, "KappaInf");
//    delete WORK;
//    delete oWork;
    return 0;
}

int improveSolution(DTDoubleArray A, DTDoubleArray B)
{
    // improve the solution
    DTMutableDoubleArray mA = A.Copy();
    DTMutableDoubleArray mB = B.Copy();
    int INFO = 3;
    int M = mA.m();
    int N = mA.n();
    int NRHS = mB.n();
    int LDA = M;
    int LDB = mB.m();
    int minDim = min(M,N);
    int IPIV[minDim];
    
    LAPACK_dgetrf(&M,&N,mA.Pointer(),&LDA,IPIV,&INFO);
    LAPACK_dgetrs(&TRANS,&N,&NRHS,mA.Pointer(),&LDA,IPIV,mB.Pointer(),&LDB,&INFO);

    int LDAF = mA.m();
    int LDX = mB.m();
    double FERR[NRHS];
    double BERR[NRHS];
    double WORK[3 * N];
    int IWORK[N];
    LAPACK_dgerfs(&TRANS, &N, &NRHS, A.Pointer(), &LDA, mA.Pointer(), &LDAF, IPIV, B.Pointer(), &LDB, mB.Pointer(), &LDX, FERR, BERR, WORK, IWORK, &INFO);
    if(INFO == 0)
    {
        DTMatlabDataFile outputFile("Output.mat");
        outputFile.Save(mB, "ximproved");
        outputFile.Save(*FERR, "Ferr");
        outputFile.Save(*BERR, "Berr");
/*
      std::cout << "ximproved is: " << std::endl;
      for(int i = 0; i < mB.m(); ++i){
	for(int j = 0; j < mB.n(); ++j)
	  {
	    std::cout << mB(i,j) << " ";
	  }
	std::cout << std::endl;
      }
      std::cout << "Ferr: " <<std::endl;
      for(int i = 0; i < NRHS; ++i)
	{
	  std::cout << FERR[i] << " ";
	}
      std::cout << std::endl;
      std::cout << "Berr: " <<std::endl;
      for(int i = 0; i < NRHS; ++i)
	{
	  std::cout << BERR[i] << " ";
	}
      std::cout << std::endl;
    }
*/
    }
//    delete [] FERR;
//    delete [] BERR;
//    delete [] WORK;
    
    return 0;
}

int main(int argc, const char *argv[])
{
    DTSetArguments(argc,argv);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    int     INFO=3;

    DTMatlabDataFile inputFile("Input.mat", DTFile::ReadOnly);
    DTDoubleArray A = inputFile.ReadDoubleArray("A");
    DTDoubleArray B = inputFile.ReadDoubleArray("b");

    if(A.n() != B.m() || A.m() == 0 || A.n() == 0 || B.m() == 0 || B.n() == 0)
    {
        std::cout << "Invalid input array dimensions." << std::endl;
        return -1;
    }

    if(calculateFactorization(A, B) != 0)
    {
        return -1;
    }

    if(solveSystem(A, B) != 0)
    {
        return -1;
    }

    if(estimateConditionNumber(A) != 0)
    {
        return -1;
    }

    if(improveSolution(A,B) != 0)
    {
        return -1;
    }
    std::cout << "program terminated."  << std::endl;
    return 0;
}  
