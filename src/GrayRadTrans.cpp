//   HyKiCT - 1D Lagrangian Radiation-Hydrodyanmics code for testing Coupling with Kinetic codes.
//   Copyright (C) 2020- Abetharan Antony
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#include "HyKiCT/GrayRadTrans.h"

GrayRadTrans::GrayRadTrans(int nx)
{
    mNx = nx;      
    mEnergyDensities.resize(mNx);
    mOpticalDepth.resize(mNx, 0);
    mgff = 1.2; //gaunt factor
    mMaterialInterface = mNx; 
}
void GrayRadTrans::setMaterialInterface(int index)
{
    mMaterialInterface = index;
}
void GrayRadTrans::storeOpacityTables(std::vector<LoadInTables> table)
{
   mOpacityTables = table; //std::move(table);
}
void GrayRadTrans::initEnergyDensity(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadEnergyDensity[i] = fixedData.radiationConstant * pow(gridData->TemperatureE[i], 4);
    }
}
void GrayRadTrans::radationPressure(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadPressure[i] = gridData->RadEnergyDensity[i] / 3; 
        gridData->TotalPressure[i] += gridData->RadPressure[i];
    }
}

void GrayRadTrans::calculateOpacities(GridData *gridData, bool calculateRoss, bool skip_search)
{
    int oi = 0;
    if(!skip_search)
    {
        mOpacityTables.at(oi).findAllOPACITYIndicies(gridData);
    }
    for(int i = 0; i < mNx; i++)
    {
        if((mMaterialInterface <= i) && (oi == 0))
        {
            oi = 1;
            if(!skip_search)
            {
                mOpacityTables.at(oi).findAllOPACITYIndicies(gridData);
            }
        }
        if(!calculateRoss)
        {
            planckianOpacity(gridData, oi, i);
            // planckianEmissionOpacity(gridData, oi, i);
        }
        else
        {
            rosslandOpacity(gridData, oi, i);
        }
    }
}

void GrayRadTrans::planckianOpacity(GridData *gridData, int table, int i)
{
    int ncols = mOpacityTables.at(table).TOPSTemperature.size();
    int elower_temp_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);
    int eupper_temp_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);

    int elower_density_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    int eupper_density_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    
    double eq11 = mOpacityTables.at(table).planckAbsorptionOpacity[elower_density_index * ncols + elower_temp_index]; //q11
    double eq12 = mOpacityTables.at(table).planckAbsorptionOpacity[elower_density_index * ncols + eupper_temp_index]; //q12
    double eq21 = mOpacityTables.at(table).planckAbsorptionOpacity[eupper_density_index * ncols + elower_temp_index]; //q21
    double eq22 = mOpacityTables.at(table).planckAbsorptionOpacity[eupper_density_index * ncols + eupper_temp_index]; //q22
    if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values) > gridData->TemperatureE[i])
    {
        gridData->RadPlanckAbsorptionOpacity[i] = 0;
    }
    else if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values) > gridData->Density[i])
    {
        gridData->RadPlanckAbsorptionOpacity[i] = 0;
    }
    else
    {
        gridData->RadPlanckAbsorptionOpacity[i] = mOpacityTables.at(table).biLinearInterpolation(std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);

    }
    gridData->RadPlanckEmissionOpacity[i] = gridData->RadPlanckAbsorptionOpacity[i]; //Set them to the same IF different opacity is present remove and add another function. 
}
void GrayRadTrans::rosslandOpacity(GridData *gridData,int table, int i)
{
    int ncols = mOpacityTables.at(table).TOPSTemperature.size();
    int elower_temp_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);
    int eupper_temp_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);

    int elower_density_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    int eupper_density_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    
    double eq11 = mOpacityTables.at(table).rosslandOpacity[elower_density_index * ncols + elower_temp_index]; //q11
    double eq12 = mOpacityTables.at(table).rosslandOpacity[elower_density_index * ncols + eupper_temp_index]; //q12
    double eq21 = mOpacityTables.at(table).rosslandOpacity[eupper_density_index * ncols + elower_temp_index]; //q21
    double eq22 = mOpacityTables.at(table).rosslandOpacity[eupper_density_index * ncols + eupper_temp_index]; //q22

    if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values) > gridData->TemperatureE[i])
    {
        gridData->RadRossAbsorptionOpacity[i] = 0;
    }
    else if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values) > gridData->Density[i])
    {
        gridData->RadRossAbsorptionOpacity[i] = 0;
    }
    else
    {   
        gridData->RadRossAbsorptionOpacity[i] = mOpacityTables.at(table).biLinearInterpolation(std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);
    }
}
void GrayRadTrans::ficksDiffusionCoefficient(GridData *gridData, FixedData const &fixedData)
{
    float larsen_limiter = fixedData.larsen_limiter;
    double diffusion, averaged_rossland_mean, rad_gradient, cell_wall_rad_energy, flux_limited_free_streaming;
    for(int i = 1; i < mNx; i++)
    {
        averaged_rossland_mean = 0.5*(gridData->Density[i - 1]* gridData->RadRossAbsorptionOpacity[i - 1] + 
                                    gridData->Density[i]* gridData->RadRossAbsorptionOpacity[i]);
        //NOTE SQUARE ROOT FLUX LIMITER Morel, J. E. 2000, J. Quant. Spectrosc. Radiat. Transfer, 65, 769. 
        rad_gradient = (gridData->RadEnergyDensity[i] - gridData->RadEnergyDensity[i - 1]) / gridData->Dx_k[i - 1]; 
        cell_wall_rad_energy = (gridData->Dx_k_1_2[i - 1] * gridData->RadEnergyDensity[i] +
                                       gridData->Dx_k_1_2[i] * gridData->RadEnergyDensity[i - 1]) / 
                                       (gridData->Dx_k_1_2[i - 1] + gridData->Dx_k_1_2[i]);
                    

        flux_limited_free_streaming = fixedData.radflxlm  * (pow(abs(rad_gradient), larsen_limiter) / pow(cell_wall_rad_energy, larsen_limiter));
        if(averaged_rossland_mean < 1.0e-50)
        {
            diffusion =  0.0; //1e-10; // This will prevent Nans in the matrix. As long as this value is << dx_k_1_2 we are safe, if you go lower density make sure to reflect that here. 
        }
        else
        {
            double larsen = 1/larsen_limiter;
            double opacity_larsen = pow(3.0 * averaged_rossland_mean, larsen_limiter);
            double bottom = pow(( opacity_larsen + flux_limited_free_streaming), larsen);
            diffusion = ((fixedData.SPEED_OF_LIGHT)) / bottom; //pow((pow(3 * averaged_rossland_mean, larsen_limiter) + flux_limited_free_streaming), 1/larsen_limiter);
        }
        gridData->RadDiffusionCoefficient[i] = diffusion; 
    }
}
void GrayRadTrans::mLeftReflectiveBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData)
{
    int i = 0;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    // mLeftBoundaryCol[0] = i ; mLeftBoundaryCol[1] = i + 1;
    std::get<0>(mLeftBoundaryValue) = 1 + r * gridData->RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]));
    std::get<1>(mLeftBoundaryValue) = -1*r * (gridData->RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
    mEnergyDensities[i] = gridData->RadEnergyDensity[i];
}
void GrayRadTrans::mRightReflectiveBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData)
{
    int i = mNx - 1;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    std::get<0>(mRightBoundaryValue) = -1*r * (gridData->RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
    std::get<1>(mRightBoundaryValue) = 1 + r * gridData->RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]));
    mEnergyDensities[i] = gridData->RadEnergyDensity[i];
}
void GrayRadTrans::mLeftVacuumSourceBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in)
{
    double D_0;
    if(gridData->RadRossAbsorptionOpacity[0] < 1.0e-50)
    {
        D_0 = 1.0e30;
    }
    else
    {
        D_0 = fixedData.SPEED_OF_LIGHT / (3.0 * gridData->Density[0]* gridData->RadRossAbsorptionOpacity[0]); 
    }
    
    double D_eff =  1 / (1 + (4 / (fixedData.SPEED_OF_LIGHT * gridData->Dx_k_1_2[0])) * D_0);
    double diff_1 = (timeData.Dt05 * gridData->RadDiffusionCoefficient[1]) / (gridData->Dx_k_1_2[0] * gridData->Dx_k[0]); 
    double diff_0 = (2*timeData.Dt05 * D_0 * D_eff) / (gridData->Dx_k_1_2[0] * gridData->Dx_k_1_2[0]); 
    int i = 0;
    double r = (timeData.Dt05 / gridData->Dx_k[i]);
    std::get<0>(mLeftBoundaryValue) = 1.0 + diff_0 + diff_1; 
    std::get<1>(mLeftBoundaryValue) =  -1.0 * diff_1;
    mEnergyDensities[i] = gridData->RadEnergyDensity[i] +  ((4.0 * diff_0) / (fixedData.SPEED_OF_LIGHT)) * f_in;

}
// void GrayRadTrans::mRightVacuumSourceBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in)
// {
//     int i = mNx - 1;
//     double D_0 = fixedData.SPEED_OF_LIGHT / (3 * gridData->Density[i]* gridData->RadRossAbsorptionOpacity[i]); 
//     double D_eff =  1 / (1 + (4 / (fixedData.SPEED_OF_LIGHT * gridData->Dx_k_1_2[i])) * D_0);
//     double diff_1 = (timeData.Dt05 * gridData->RadDiffusionCoefficient[mNx - 1]) / (gridData->Dx_k_1_2[i] * gridData->Dx_k[i]); 
    // double diff_0 = (2*timeData.Dt05 * D_0 * D_eff) / (gridData->Dx_k_1_2[i] * gridData->Dx_k_1_2[i]); 
    // double r = (timeData.Dt05 / gridData->Dx_k[i]);
    // std::get<0>(mRightBoundaryValue) = 1 + diff_0 + diff_1 + timeData.Dt05; 
    // std::get<1>(mRightBoundaryValue) =  -1 * diff_1;
    // mEnergyDensities[i] = gridData->RadEnergyDensity[i] +  ((4 * diff_0) / (fixedData.SPEED_OF_LIGHT)) * f_in;
// }
void GrayRadTrans::ficksDiffusion(GridData *gridData, FixedData const &fixedData, TimeStep const &timeData, MatrixSolver &matrixSolver)
{   //Method from Bowers And Wilson p427 .. modified for my needs 
    //See lab book for currently implemented method. 
    //Matrix Equations being solved is:
    // (1 + r(Dhat_{k+1} + Dhat_{k}))U^{n+1}_{k+1/2} - rDhat_{k+1}U^{n+1}_{k+3/2} - rDhat_{k}U^{n+1}_{k-1/2} = U^{''}_{k+1/2}
    // r = dt_{n+1/2} / dr{k+1/2}
    // Dhat_{k+1} = D_r / dr{k+1}
    // D_r_k = (dr * SPEED_OF_LIGHT/3) / (rho_k * Kappa_Ross * dr_k + C_L mag(u_{k+1/2} - u_{k-1/2}) /0.5(u_{k+1/2} - u_{k-1/2}))
    // GrayRadTrans::createPETScObjects();
    std::vector<std::tuple<double, double, double>> diff_coef;
    bool direchletRight = false, direchletLeft = false;
    double val_0,val_1,val_2;
    diff_coef.resize(mNx);
    for(int i = 1; i < mNx - 1; i++)
    {
        double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
        val_0 = -1*r * (gridData->RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
        val_1 = 1 + r * gridData->RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]))
                    + r * gridData->RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]));
        val_2 = -1*r * (gridData->RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
        diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
        mEnergyDensities[i] = gridData->RadEnergyDensity[i];
    }
    if(fixedData.RightRadBoundaryCondition == "r")
    {
        mRightReflectiveBoundaryCondition(gridData, fixedData, timeData);
    }
    else
    {
        direchletRight = true;
        mEnergyDensities[mNx - 1] = fixedData.RightDirichletEnergyDensity;
    }
    if(fixedData.LeftRadBoundaryCondition == "r")
    {
        mLeftReflectiveBoundaryCondition(gridData, fixedData, timeData);
    }
    else if((fixedData.LeftRadBoundaryCondition == "v") || (fixedData.LeftRadBoundaryCondition == "m"))
    {
        double f_in = 0.0;
        if(fixedData.LeftRadBoundaryCondition  == "m")
        {
            f_in =  (fixedData.SPEED_OF_LIGHT / 4.0) * fixedData.radiationConstant * pow(fixedData.LeftSourceTemperature, 4);
        }

        mLeftVacuumSourceBoundaryCondition(gridData, fixedData, timeData, f_in);
    }
    else
    {
        direchletLeft = true;
        mEnergyDensities[0] = fixedData.LeftDirichletEnergyDensity;
    }
    // mLeftBoundaryCondition(gridData, fixedData, timeData);
    //     mRightBoundaryCondition(gridData, fixedData, timeData);

    matrixSolver.fillTriMatrix(diff_coef, mEnergyDensities, mLeftBoundaryValue, mRightBoundaryValue,
                                 direchletLeft,direchletRight);
    matrixSolver.solveMatrix();
    matrixSolver.returnSolution(gridData->RadFicks);
}

void GrayRadTrans::updateRadiationTemperature(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadTemperature[i] = pow(((fixedData.SPEED_OF_LIGHT /  (4 * fixedData.STEFAN_BOLTZMANN_CONSTANT)) * gridData->RadEnergyDensity[i]), 0.25);
    }
}

void GrayRadTrans::calcFreeFreeAbsorption(GridData *gridDataT,FixedData const&fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->RadFFAbsorb[i] = fixedData.SPEED_OF_LIGHT * gridDataT->RadPlanckAbsorptionOpacity[i] *gridDataT->RadEnergyDensity[i];
        assert(gridDataT->RadFFAbsorb[i]>=0);
    }
}
void GrayRadTrans::calcFreeFreeEmission(GridData *gridDataT,FixedData const&fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->RadFFEmission[i] = gridDataT->RadPlanckEmissionOpacity[i] *
                                        4.0 * fixedData.STEFAN_BOLTZMANN_CONSTANT * pow(gridDataT->TemperatureE[i], 4);
    }
}
void GrayRadTrans::updateRadiationEnergyDensity(GridData *gridDataT,FixedData const&fixedData, TimeStep const &timeData)
{      
    // Update for emission
    // Emission(J) - Absorption(A) WKg-1
    // U'' = U' + rho*(4 *sigma* kpE(Te) * Te^4 - ckpA(Te, Tr) *Ur)
    calcFreeFreeAbsorption(gridDataT, fixedData);
    calcFreeFreeEmission(gridDataT, fixedData);

    for(int i = 0; i < mNx; i++)
    {
        gridDataT->RadEnergyDensity[i] += gridDataT->Density[i] * timeData.Dt05 * (gridDataT->RadFFEmission[i] - gridDataT->RadFFAbsorb[i]);
    }
}
void GrayRadTrans::updateRadiationEnergyDensity(GridData *gridDataT)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->RadEnergyDensity[i] = gridDataT->RadFicks[i];
    }
}

void GrayRadTrans::updateRadiationQuantities(GridData *gridDataT, FixedData const &fixedData, TimeStep const &timeData, MatrixSolver &matrixSolver, int updateType)
{
    if(updateType == 1)
    {   /*Update for emission
        Update Chain
        Calculate Rad Coefficient
        Calculate new heat Capacity which accounts for the radiation emission. 
        Calculate Tr'' 
        Update U''
        U'' = U' + J - A*/
        calculateOpacities(gridDataT, false);
        updateRadiationEnergyDensity(gridDataT, fixedData, timeData);
        updateRadiationTemperature(gridDataT, fixedData);
    }
    else if(updateType == 2)
    {   //Update for Ficks Diffusion
        //Update Chain
        //Calculate opacities
        //Calcualte Diffusion coefficient
        //Calculate the diffusion to find U^n+1
        //U^(n+1) = U_ficks
        calculateOpacities(gridDataT, true, true);
        // RosslandFreeFree(gridDataT, fixedData);
        double sum_check = std::accumulate(gridDataT->RadRossAbsorptionOpacity.begin(), gridDataT->RadRossAbsorptionOpacity.end(), 0.0);
        if(sum_check > 1e-50)
        {
            ficksDiffusionCoefficient(gridDataT, fixedData);
            ficksDiffusion(gridDataT, fixedData, timeData, matrixSolver);
            updateRadiationEnergyDensity(gridDataT);
            // calcFreeFreeEmission(gridDataT, fixedData);
            // calcFreeFreeAbsorption(gridDataT, fixedData);
            updateRadiationTemperature(gridDataT, fixedData);
        }
    }
}
void GrayRadTrans::setPlanckianOpacity(GridData *gridData, double value)
{
    std::fill(gridData->RadPlanckEmissionOpacity.begin(), gridData->RadPlanckEmissionOpacity.end(),value); 
    std::fill(gridData->RadPlanckAbsorptionOpacity.begin(), gridData->RadPlanckAbsorptionOpacity.end(),value); 
}
void GrayRadTrans::setRossOpacity(GridData *gridData, double value)
{
    std::fill(gridData->RadRossAbsorptionOpacity.begin(), gridData->RadRossAbsorptionOpacity.end(), value);
}
void GrayRadTrans::setRadEnergyDensity(GridData *gridData, double value)
{
    std::fill(gridData->RadEnergyDensity.begin(), gridData->RadEnergyDensity.end(), value); 
}
void GrayRadTrans::setRadEnergyDensity(GridData *gridData, std::vector<double> values)
{
    for(int i = 0; i < values.size(); i++)
    {
        gridData->RadEnergyDensity[i] = values[i];
    }
}