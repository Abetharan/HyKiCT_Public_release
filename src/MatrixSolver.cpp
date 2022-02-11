#include "HyKiCT/MatrixSolver.h"
MatrixSolver::MatrixSolver(int Nx)
{
    mNx = Nx;
    mn = Nx;
}
MatrixSolver::~MatrixSolver()
{
    // cleanUpPETScObjects();
} 
void MatrixSolver::cleanUpPETScObjects()
{
    VecDestroy(&mb);
    VecDestroy(&mSol); MatDestroy(&mA);   
    KSPDestroy(&mKSP);
}
void MatrixSolver::createPETScObjects()
{
    VecCreate(PETSC_COMM_WORLD, &mSol);
    PetscObjectSetName((PetscObject) mSol, "Solution");
    VecSetSizes(mSol, PETSC_DECIDE, mn);
    VecSetFromOptions(mSol);
    
    VecCreate(PETSC_COMM_WORLD, &mb);
    PetscObjectSetName((PetscObject) mb, "RHS");
    VecSetSizes(mb, PETSC_DECIDE, mn);
    VecSetFromOptions(mb);
    //VecDuplicate(msol, &mb);
    // PetscObjectSetName((PetscObject) mb, "RHS");
    //Create Petsc Mat
    // std::cout << "Creating MATRICES" << '\n';
    MatCreate(PETSC_COMM_WORLD, &mA);
    MatSetSizes(mA, PETSC_DECIDE, PETSC_DECIDE, mn, mn); //size set to n = mNx
    PetscObjectSetName((PetscObject) mA, "Matrix");
    MatSetFromOptions(mA);
    MatSetUp(mA);
    KSPCreate(PETSC_COMM_WORLD, &mKSP);
    KSPSetType(mKSP, "bcgs");
    KSPGetPC(mKSP, &mPC);
    
    PCSetType(mPC, PCJACOBI);
    PCFactorSetShiftType(mPC, MAT_SHIFT_NONZERO);
}
void MatrixSolver::fillTriMatrix(std::vector<std::tuple<double, double, double>> D, std::vector<double> S,
                         std::tuple<double, double> leftBC, std::tuple<double, double> rightBC, bool leftDirechlet,
                         bool rightDirechlet)
{
    for(i = 1; i < mn - 1; i++)
    {
        mCol[0] = i - 1; mCol[1] = i; mCol[2] = i + 1;

        mValue[0] = std::get<0>(D[i]);
        mValue[1] = std::get<1>(D[i]);
        mValue[2] = std::get<2>(D[i]);
        mSource = (PetscReal)(S[i]);
        MatSetValues(mA, 1, &i, 3, mCol, mValue, INSERT_VALUES); 
        VecSetValues(mb, 1, &i, &mSource, INSERT_VALUES);
    }
    if(!leftDirechlet)
    {
        i = 0;
        mLeftBoundaryCol[0] = i ; mLeftBoundaryCol[1] = i + 1;
        mLeftBoundaryValue[0] = std::get<0>(leftBC);
        mLeftBoundaryValue[1] = std::get<1>(leftBC);
        mSource = (PetscReal)(S[i]);
        MatSetValues(mA, 1, &i, 2, mLeftBoundaryCol, mLeftBoundaryValue, INSERT_VALUES);
        VecSetValues(mb, 1, &i, &mSource, INSERT_VALUES);
    }
    else
    {
        i = 0;  mStartCol = 0;
        mStartValue = 1;
        mSource = (PetscReal)(S[i]);
        MatSetValues(mA, 1, &i, 1, &mStartCol, &mStartValue, INSERT_VALUES);
        VecSetValues(mb, 1, &i, &mSource, INSERT_VALUES);
    }
    if(!rightDirechlet)
    {
        i = mNx - 1;
        mRightBoundaryCol[0] = i - 1; mRightBoundaryCol[1] = i;
        mRightBoundaryValue[0] = std::get<0>(rightBC);
        mRightBoundaryValue[1] = std::get<1>(rightBC);
        mSource = (PetscReal)(S[i]);
        MatSetValues(mA, 1, &i, 2, mRightBoundaryCol, mRightBoundaryValue, INSERT_VALUES);
        VecSetValues(mb, 1, &i, &mSource, INSERT_VALUES);

    }
    else
    {
        i = mn - 1; mLastCol= mn - 1;
        mLastValue = 1;
        mSource = (PetscReal)(S[i]);
        MatSetValues(mA, 1, &i, 1, &mLastCol, &mLastValue, INSERT_VALUES);
        VecSetValues(mb, 1, &i, &mSource, INSERT_VALUES);
    }
    VecAssemblyBegin(mb); 
    MatAssemblyBegin(mA, MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(mb);
    MatAssemblyEnd(mA, MAT_FINAL_ASSEMBLY);
    // std::cout << "OUTPUTTING VEC AND MATRIX" << '\n';
    // MatView(mA, PETSC_VIEWER_STDOUT_WORLD);
    // VecView(mb, PETSC_VIEWER_STDOUT_WORLD);
}
void MatrixSolver::solveMatrix()
{
    
    VecNorm(mb, NORM_2, &mNorm);
    KSPSetTolerances(mKSP, 1e-20, 1e-15 * mNorm, 1e7, 1000);
    // KSPSetFromOptions(mKSP);
    
    KSPSetOperators(mKSP, mA, mA);
    KSPSolve(mKSP, mb, mSol);
    // KSPView(mKSP, PETSC_VIEWER_STDOUT_WORLD);

    //VecAXPY(msol, -1.0, mu);
    KSPGetConvergedReason(mKSP, &mKSPReason);
    KSPGetIterationNumber(mKSP, &mits);
    // PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g, Iterations %D\n", (double)mNorm, mits);
    //VecView(msol, PETSC_VIEWER_STDERR_WORLD);
    if(mKSPReason == KSP_DIVERGED_INDEFINITE_PC)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\nDivergence because of indefinite preconditioner;\n");
        PetscPrintf(PETSC_COMM_WORLD, "Run the executable again but with '-pc facto shift type POSITIVE_DEFININTE opition.\n");
    }
    else if(mKSPReason < 0)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\nOther kind of divergene: this should not happen.\n");
    }
    else
    {
        KSPGetIterationNumber(mKSP, &mits); 
    }
}
void MatrixSolver::returnSolution(std::vector<double> &fillVector)
{
    PetscScalar const *value_1;
    i = 0;
    VecGetArrayRead(mSol, &value_1);
    for(i = 0; i < mNx; i++)
    {
        // VecGetValues(mSol, 1, &i, &value_1);
        fillVector[i] = value_1[i];
    }
    VecRestoreArrayRead(mSol, &value_1);
}
