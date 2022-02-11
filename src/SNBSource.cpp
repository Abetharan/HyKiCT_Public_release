#include "HyKiCT/SNBSource.h"

SNBSource::SNBSource(int nx)
{
    mNx = nx;
}
void SNBSource::InitSNB(int ng, int id, double maxE)
{
    mNg = ng;
    mGroupID = id;
    mDqNL.resize(mNx - 1);
    mEField.resize(mNx - 1);
    mSNBGridWall.resize(mNg + 1);
    mSNBGridCentered.resize(mNg);
    mSpitzer.resize(mNx + 1);
    mLambdaStarWall.resize(mNx - 1);
    mLambdaStarCentre.resize(mNx);
    mLambdaE.resize(mNx - 1);
    mSNBBetaCentre.resize(mNx);
    mSNBDbetaCentre.resize(mNx);
    mSNBBetaWall.resize(mNx - 1);
    mSNBDbetaWall.resize(mNx - 1);
    mSNBH.resize(mNx);
    mSNBU.resize(mNx - 1);
    for(int i = 0; i < mNx; i++)
    {
        mLambdaStarCentre[i].resize(mNg);
        if (i < mNx - 1)
        {
            mLambdaE[i].resize(mNg);
            mLambdaStarWall[i].resize(mNg);
            mSNBBetaWall[i].resize(mNg);
            mSNBDbetaWall[i].resize(mNg);
            mSNBU[i].resize(mNg); 
        }
        mSNBBetaCentre[i].resize(mNg);
        mSNBDbetaCentre[i].resize(mNg);
        mSNBH[i].resize(mNg);
    }
    if (mNg == 1) 
    {
        mSNBGridWall.push_back(0);
    }

    double delta = (maxE - 0) / (mNg);
    for(int g = 0; g < mNg; ++g)
    {
        mSNBGridWall[g] = (0 + delta * g);
        if(g > 0)
        {
            mSNBGridCentered[g - 1] = ((mSNBGridWall[g] + mSNBGridWall[g - 1]) / 2);
        }
    }
    mSNBGridWall.back() = maxE; 
    mSNBGridCentered.back() = (mSNBGridWall.back() + *(mSNBGridWall.end() - 2))/ 2;
}

void SNBSource::mCalculateSNBBeta(std::vector<double> const &TemperatureE)
{ //relies on interpolation
    for(int i = 0; i < mNx; i++)
    {
        int g = mGroupID;
        {
            mSNBBetaCentre[i][g] = (mSNBGridCentered[g]/TemperatureE[2*i]);
            mSNBDbetaCentre[i][g] = (mSNBGridWall[g + 1] - mSNBGridWall[g]) / TemperatureE[2*i];
            if(i < mNx - 1)
            {
                mSNBBetaWall[i][g] = (mSNBGridCentered[g]/TemperatureE[2*i + 1]);
                mSNBDbetaWall[i][g] = (mSNBGridWall[g + 1] - mSNBGridWall[g]) / TemperatureE[2*i + 1];
            }
        }
    }
}
void SNBSource::mCalculateSNBLambdaStar(std::vector<double> const &TemperatureE, std::vector<double> const &NumberDensityE,
                             std::vector<double> const &Zbar, std::vector<double> const &CoulombLog,
                             FixedData const &fixedData)
{
    // double gamma_ee = pow(fixedData.ELECTRON_CHARGE, 2) / (4 * M_PI * fixedData.VACUUM_PERMITTIVITY * 
    //                     fixedData.ELECTRON_MASS);
    double gamma_ee = (fixedData.ELECTRON_CHARGE * fixedData.ELECTRON_CHARGE) / (4 * M_PI * fixedData.VACUUM_PERMITTIVITY * 
                        fixedData.ELECTRON_MASS);
    double Y = 4 * M_PI * gamma_ee * gamma_ee;
    int j = 0; int k = 0;
    double vt_centre, tau_ei_centre, zeta_centre;
    double vt_wall, tau_ei_wall, zeta_wall;
    for(int i = 0; i < mNx; i++)
    {
        //Note centre quantities are used in the constant factor part of the equation 
        vt_centre = pow(((2*fixedData.BOLTZMANN_CONSTANT * TemperatureE[2*i]) / 
                            fixedData.ELECTRON_MASS), 0.5);
        // tau_ei_centre = pow(vt_centre, 3) / (Y * Zbar[2*i] * NumberDensityE[2*i] 
        //                 * CoulombLog[2 * i]);
        tau_ei_centre = (vt_centre*vt_centre*vt_centre) / (Y * Zbar[2*i] * NumberDensityE[2*i] 
                        * CoulombLog[2 * i]);
        if(mSeperated) //Implements the suggested method in brodrick
        {
            zeta_centre = Zbar[2*i] / 2 ;
        }
        else //Implements the corrected averaged value as suggested by brodrick
        {
            zeta_centre = pow(Zbar[2*i] / 2, 0.5) *InvSqrt((Zbar[2 * i] + 4.2) / (Zbar[2 * i] + 0.24)); 

        }
        double lambda_ei_centre = vt_centre * tau_ei_centre;
        double lambda_e_centre = zeta_centre * lambda_ei_centre;

        if(i < mNx - 1)
        {   //Note wall quantities are used in the diffusion part of the equation 
            vt_wall= pow(((2*fixedData.BOLTZMANN_CONSTANT * TemperatureE[2 * i + 1]) / fixedData.ELECTRON_MASS), 0.5);
            // tau_ei_wall= pow(vt_wall, 3) / (Y * Zbar[2 * i + 1] * NumberDensityE[2 * i + 1] * CoulombLog[2 * i + 1]);
            tau_ei_wall= (vt_wall*vt_wall*vt_wall) / (Y * Zbar[2 * i + 1] * NumberDensityE[2 * i + 1] * CoulombLog[2 * i + 1]);
            if(mSeperated) //Implements the suggested method in brodrick
            {
                zeta_wall = (Zbar[2 * i + 1] + 0.24) / (Zbar[2 * i + 1] + 4.2);
            }
            else //Implements the corrected averaged value as suggested by brodrick
            {
                zeta_wall= pow(Zbar[2*i + 1] / 2, 0.5) * InvSqrt((Zbar[2 * i + 1] + 4.2) / (Zbar[2 * i + 1] + 0.24));

            }
        }
        double lambda_ei_wall = vt_wall* tau_ei_wall;
        double lambda_e_wall= zeta_wall* lambda_ei_wall;

        // for(int g = 0; g<mNg; g++)
        int g = mGroupID;
        {
            mLambdaStarCentre[i][g] = (mSNBBetaCentre[i][g] * mSNBBetaCentre[i][g]) * lambda_e_centre;
            if(i < mNx - 1)
            {
                mLambdaStarWall[i][g] = (mSNBBetaWall[i][g]*mSNBBetaWall[i][g]) * lambda_e_wall; 
            }
        }
    }
}

void SNBSource::mCalculateSpitzerElectricField(std::vector<double> const &cell_centre_grid, std::vector<double> const &TemperatureE, 
                                    std::vector<double>  const &NumberDensityE, std::vector<double> const &Zbar, FixedData const &fixedData)
{
    for(int i = 1; i < mNx; i++)
    {
        double dx = cell_centre_grid[i] - cell_centre_grid[i - 1];
        double gamma =  1 + ((3 * (Zbar[i] + 0.477)) / (2*(Zbar[i] + 2.15)));
        mEField[i - 1] = -1 * (fixedData.BOLTZMANN_CONSTANT/fixedData.ELECTRON_CHARGE) * TemperatureE[2*i] * 
                     (((NumberDensityE[2*i] - NumberDensityE[2*i - 2]) / (NumberDensityE[2*i - 1] * dx)) +
                    gamma * ((TemperatureE[2*i] - TemperatureE[2*i - 2]) /(TemperatureE[2*i - 1] * dx)));
    }
}
void SNBSource::mCalculateSNBLambdaE(FixedData const &fixedData)
{
    double inverse_lambda;
    for(int i = 0; i < mNx - 1; i++)
    {
        int g = mGroupID;
        {  
            inverse_lambda =  1 / (mLambdaStarWall[i][g]) + 
                abs((fixedData.ELECTRON_CHARGE * mEField[i]) / (fixedData.BOLTZMANN_CONSTANT * mSNBGridCentered[g]));
            mLambdaE[i][g] = 1/inverse_lambda;
        }
    }
}
void SNBSource::mCalculateSNBU(std::vector<double> const &SpitzerHarmHeatFlow)
{
   for(int i = 1; i < mNx; i++)
    {
        int g = mGroupID;
        // for(int g = 0; g < mNg; g++)
        {                       //Recall that Our Spitzer harm heat flow is opposite sign
            // mSNBU[i - 1][g] = (-1*SpitzerHarmHeatFlow[i] / 24) * mSNBDbetaWall[i - 1][g] * pow(mSNBBetaWall[i - 1][g], 4) * exp(-1*mSNBBetaWall[i - 1][g]);
            mSNBU[i - 1][g] = (-1*SpitzerHarmHeatFlow[i] / 24) * mSNBDbetaWall[i - 1][g] * (mSNBBetaWall[i - 1][g]*mSNBBetaWall[i - 1][g]*mSNBBetaWall[i - 1][g]*mSNBBetaWall[i - 1][g]) * exp(-1*mSNBBetaWall[i - 1][g]);
        }
    }
}
void SNBSource::mCalculateSNBH(GridData const *gridData, std::vector<double> const &Zbar)
{  
    double dx_diff, dx_wall_back, dx_wall_forward, source;
    // std::vector<std::tuple<double, double, double>> diff_coef;
    double val_0,val_1,val_2;
    // std::vector<double> source_vec;
    // source_vec.resize(mNx);
    // diff_coef.resize(mNx);
    double a[mNx];
    double b[mNx];
    double c[mNx]; 
    double d[mNx];
    // std::vector<double> aa(mNx -1 ),bb(mNx- 1),cc(mNx - 1),dd(mNx -1);
    // for(int g = 0; g < mNg; g++)
    int g = mGroupID;
    {
        for(int i = 1; i < mNx - 1; i++)
        {
            // dx_diff = cell_wall_grid[i] - cell_wall_grid[i - 1];
            // dx_wall_back = cell_centre_grid[i] - cell_centre_grid[i - 1];
            // dx_wall_forward = cell_centre_grid[i + 1] - cell_centre_grid[i];
            dx_diff = gridData->Dx_k_1_2[i];
            dx_wall_back = gridData->Dx_k[i - 1];
            dx_wall_forward = gridData->Dx_k[i];
            val_0 = -1 * mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff);
            val_1 = mLambdaE[i][g] / (3*dx_wall_forward * dx_diff) + 
                        mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff) +
                        mR / (mLambdaStarCentre[i][g]) ; //*Zbar[i]
            val_2 = -1*mLambdaE[i][g] / (3*dx_wall_forward * dx_diff);
            // source_vec[i] = -1 * ((mSNBU[i][g] - mSNBU[i - 1][g] ) / dx_diff);
            // diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
            a[i] = val_0;
            b[i] = val_1;
            c[i] = val_2;
            d[i] = -1 * ((mSNBU[i][g] - mSNBU[i - 1][g] ) / dx_diff);
            // aa[i - 1] = val_0;
            // bb[i] = val_1;
            // cc[i] = val_2;
            // dd[i] = -1 * ((mSNBU[i][g] - mSNBU[i - 1][g] ) / dx_diff);
        }
        
        int i = 0;
        // dx_wall_forward = cell_wall_grid[1];
        dx_wall_forward = gridData->Dx_k[i];
        dx_diff = 2 * (gridData->CellCenteredCoord[0] - gridData->CellWallCoord[0] );
        // std::get<0>(mLeftBoundaryValue) = mLambdaE[i][g] / (3*dx_wall_forward * dx_diff) + 
        //                 mR / (mLambdaStarCentre[i][g] * Zbar[i]);
        // std::get<1>(mLeftBoundaryValue) = -1*mLambdaE[i][g] / (3*dx_wall_forward * dx_diff); 
        // source_vec[i] = -1 * (mSNBU[i][g]/ dx_diff);
        b[i] = mLambdaE[i][g] / (3*dx_wall_forward * dx_diff) + 
                        mR / (mLambdaStarCentre[i][g]); // * Zbar[i]
        c[i] =-1*mLambdaE[i][g] / (3*dx_wall_forward * dx_diff);
        d[i] = -1 * (mSNBU[i][g]/ dx_diff); 

        // bb[i] = mLambdaE[i][g] / (3*dx_wall_forward * dx_diff) + 
        //                 mR / (mLambdaStarCentre[i][g] * Zbar[i]);
        // cc[i] = -1*mLambdaE[i][g] / (3*dx_wall_forward * dx_diff);
        // dd[i] = -1 * (mSNBU[i][g]/ dx_diff);

        i = mNx - 1;
        // dx_wall_back = cell_wall_grid[mNx] - cell_wall_grid[mNx - 1];
        dx_wall_back = gridData->Dx_k[i - 1];
        dx_diff = 2 * (gridData->CellWallCoord[mNx] - gridData->CellCenteredCoord[mNx - 1]);
        // std::get<0>(mRightBoundaryValue) = -1*mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff); 
        // std::get<1>(mRightBoundaryValue) = mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff) + 
                        // mR / (mLambdaStarCentre[i][g] * Zbar[i]);
        // source_vec[i] = 1 * (mSNBU[i - 1][g]/ dx_diff);
        
        a[i] = -1*mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff);  
        b[i] = mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff) + 
                        mR / (mLambdaStarCentre[i][g] ); //* Zbar[i]
        d[i] = 1 * (mSNBU[i - 1][g]/ dx_diff); 
        // bb[i] = mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff) + 
        //                 mR / (mLambdaStarCentre[i][g] * Zbar[i]); 
        // aa[i] =-1*mLambdaE[i - 1][g] / (3*dx_wall_back * dx_diff);   
        // dd[i]=1 * (mSNBU[i - 1][g]/ dx_diff); 


        // matrixSolver.fillTriMatrix(diff_coef, source_vec, mLeftBoundaryValue, mRightBoundaryValue, false, false);
        // matrixSolver.solveMatrix();
        // std::vector<double> fillVector;
        // fillVector.resize(mNx);
        // matrixSolver.returnSolution(fillVector);
        
        mSolveTridiagonal(a, b, c, d, mNx);
        // //Extract value from PETSC vector and store into std::vector 
        // PetscScalar value_1;
        for(int j = 0; j < mNx; j++)
        {
            mSNBH[j][g] = d[j];//fillVector[j];
            
        }
    }
}
void SNBSource::mCalculateSNBDqNL(std::vector<double> &SNBDqNL, std::vector<double> const &cell_centre_grid)
{
    int g = mGroupID;
    for(int i = 0; i < mNx - 1; i++)
    {
        double heat_flow = 0;
        // for(int g = 0; i < mNg; g++)
        // {
            double dx_wall = cell_centre_grid[i + 1] - cell_centre_grid[i]; 
            double grad_delta_f0 = ((mSNBH[i + 1][g] - mSNBH[i][g])/dx_wall);
            heat_flow = mLambdaE[i][g] * grad_delta_f0;
        // }
        SNBDqNL[i] += 0.333333333333333333333 * heat_flow;
    }
}
void SNBSource::snbHeatFlowCorrection(GridData const *gridData, FixedData const &fixedData)
{
    mCalculateSNBBeta(gridData->InterTemperatureE);
    mCalculateSNBLambdaStar(gridData->InterTemperatureE, gridData->InterNumberDensityE, gridData-> InterZbar,
                            gridData->InterCoulombLog, fixedData);
    mCalculateSpitzerElectricField(gridData->CellCenteredCoord, gridData->InterTemperatureE, 
                                        gridData->InterNumberDensityE, gridData-> Zbar, fixedData);
    mCalculateSNBLambdaE(fixedData);
    mCalculateSNBU(gridData->HeatFlowE);
    // mCalculateSNBH(gridData->Dx_k, gridData->Dx_k_1_2, gridData->Zbar, matrixSolver);
    mCalculateSNBH(gridData, gridData->Zbar);
}
void SNBSource::SNBGather(GridData *gridData)
{
    mCalculateSNBDqNL(gridData->SNBDqNL, gridData->CellCenteredCoord);
}
void SNBSource::SNBGatherH(GridData *gridData)
{
    if(mGroupID == 0)
    {
        gridData->SNBH.resize(mNx);
    }
    for(int i = 0; i < mNx; i++)
    {
        if(mGroupID == 0)
        {
            gridData->SNBH[i].resize(mNg);
        }
        gridData->SNBH[i][mGroupID] = mSNBH[i][mGroupID];
    }
}
void SNBSource::SNBCorrect(GridData *gridData)
{
    for(int i = 1; i < mNx; i++)
    {
        gridData->HeatFlowE[i] += gridData->SNBDqNL[i-1];
    }
}

void SNBSource::mSolveTridiagonal(double* a, double* b, double* c, double* d, int n) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}
