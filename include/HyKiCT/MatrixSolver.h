#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <tuple>
#include "petscksp.h" 
#include "petscmat.h"
#include "petscvec.h"
#include "GridData.h"
#include "FixedData.h"
#include "Switches.h"
#include "TimeStep.h"

class MatrixSolver
{
    private:
        Vec mb,mSol;
        Mat mA;
        KSP mKSP;
        PC mPC;
        PetscReal mNorm;
        // PetscInt i,mn, mLastCol[2], mStartCol[2], mCol[3], mits;
        PetscInt i,mn, mLastCol, mStartCol, mCol[3],mLeftBoundaryCol[2], mRightBoundaryCol[2], mits;
        KSPConvergedReason mKSPReason;
        
        //PetscInitialize()
        PetscScalar mLastValue, mStartValue,mValue[3],mLeftBoundaryValue[2], mRightBoundaryValue[2], mGradU, mSource;
        int mNx;
    public:
        void createPETScObjects();
        void cleanUpPETScObjects();
        MatrixSolver(int Nx);
        ~MatrixSolver();
        void fillTriMatrix(std::vector<std::tuple<double, double, double>> D, std::vector<double> S,
                         std::tuple<double, double> leftBC,std::tuple<double, double> rightBC,
                        bool leftDirechlet,bool rightDirechlet);
        void ThomasTridiagonalolver(double* a, double* b, double* c, double* d, int n);
        void solveMatrix();
        void returnSolution(std::vector<double> &fillVector);
};
#endif