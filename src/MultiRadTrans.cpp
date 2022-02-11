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
#include "HyKiCT/MultiRadTrans.h"

MultiRadTrans::MultiRadTrans(int nx)
{
    mNx = nx;      
    mEnergyDensities.resize(mNx);
    mOpticalDepth.resize(mNx, 0);
    mgff = 1.2; //gaunt factor
    mMaterialInterface = mNx; 
    mLocal.RadPlanckAbsorptionOpacity.resize(mNx);
    mLocal.RadPlanckEmissionOpacity.resize(mNx);
    mLocal.RadRossAbsorptionOpacity.resize(mNx);
    mLocal.RadEnergyDensity.resize(mNx, 1e-15);
    mLocal.RadFFEmission.resize(mNx);
    mLocal.RadFFAbsorb.resize(mNx);
    mLocal.RadDiffusionCoefficient.resize(mNx + 1);
    mPlanck.resize(mNx);
    mPlanckIntegral.resize(mNx);
    mLeakAreaVector.resize(mNx);
    mLeakSourceCoefficient.resize(mNx);
}
void MultiRadTrans::calculateRadPressure(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->TotalPressure[i] += gridData->RadTemperature[i]*gridData->RadTemperature[i]*gridData->RadTemperature[i]*gridData->RadTemperature[i] *
                                     (4/3) * (fixedData.radiationConstant * fixedData.SPEED_OF_LIGHT * 0.25);
    }
}
double MultiRadTrans::exp2(double y) {
 
    double EXP_A = (1048576/0.69314718055994530942);
    double EXP_C = 60801;
    return EXP_A*(y) + (1072693248 - EXP_C);
}
double MultiRadTrans::mInterpolateArea(double time, double lower_time, double upper_time, double lower_power, double upper_power)
{
    double weight = (time - lower_time) / (upper_time - lower_time);
    double interIntensity = lower_power * (1 - weight) + upper_power * weight; 
    return interIntensity;
}
void MultiRadTrans::mFindLeakArea(double search_time, FixedData const &fixedData)
{
    double multiplier1, multiplier2;
    // double area1 = 8.6e-4; //assuming an inital radius of 1.5E-3;
    // double area2= 2.1e-4;
    double area1 = fixedData.leakArea1;//3.6e-4; //assuming an inital radius of 1.5E-3;
    double area2 = fixedData.leakArea2;//0.5e-4;
    // std::vector<double> fixedData.leakAreaMultiplier1 = {1.0, 0.999, 0.999, 0.8, 0.523, 0.23, 0.02, 0.01};
    // std::vector<double> fixedData.leakAreaMultiplier2 = {1.0, 1.0,1.0,1.0,1.0, 0.95, 0.55, 0.15};
    // std::vector<double> fixedData.leakTimes = {0, 14.5e-9, 17.6e-9, 19.0e-9, 20.0e-9, 21.0e-9, 22.0e-9, 23.0e-9};
    auto time_upper_bound_search = std::lower_bound(fixedData.leakTimes.begin(), fixedData.leakTimes.end(), search_time );
    auto time_lower_bound_search = time_upper_bound_search - 1;
    int time_lower_index = time_lower_bound_search - fixedData.leakTimes.begin();
    int time_upper_index = time_upper_bound_search - fixedData.leakTimes.begin();
    if(time_lower_index < 0)
    {
        time_lower_index = 0;
    }
    if((time_lower_index == 0) && (time_upper_index == 0))
    {
        multiplier1 = fixedData.leakAreaMultiplier1[0];
        multiplier2 = fixedData.leakAreaMultiplier2[0];
    }
    else
    {
        multiplier1 = mInterpolateArea(search_time, fixedData.leakTimes[time_lower_index], fixedData.leakTimes[time_upper_index],
                                    fixedData.leakAreaMultiplier1[time_lower_index], fixedData.leakAreaMultiplier1[time_upper_index]);
        multiplier2 = mInterpolateArea(search_time, fixedData.leakTimes[time_lower_index], fixedData.leakTimes[time_upper_index],
                                    fixedData.leakAreaMultiplier2[time_lower_index], fixedData.leakAreaMultiplier2[time_upper_index]);
    }
    mLeakArea =  area1*multiplier1 + area2*multiplier2;
}
void MultiRadTrans::setLeakArea(FixedData const &fixedData, double search_time)
{   
    mFindLeakArea(search_time, fixedData);
    int start_index, end_index;
    if(fixedData.leakMaterial == 0)
    {
        start_index = 0;
        end_index = mMaterialInterface;
    }
    else if(fixedData.leakMaterial == 1)
    {
        start_index = mMaterialInterface;
        end_index = mNx;
    }
    else
    {
        start_index = 0;
        end_index = mNx;
    }
    for(int i = start_index; i < end_index; i++ )
    {
        mLeakAreaVector[i] = mLeakArea;
    }
}

void MultiRadTrans::setMaterialInterface(int index)
{
    mMaterialInterface = index;
}
void MultiRadTrans::storeOpacityTables(std::vector<LoadInTables> table)
{
   mOpacityTables = table;//std::move(table);
}
void MultiRadTrans::setEnergyGroup(int group)
{
   mEnergyGroup = group;
}
void MultiRadTrans::initEnergyDensity(GridData const *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        mLocal.RadEnergyDensity[i] = fixedData.radiationConstant * pow(gridData->TemperatureE[i], 4.0);
    }
}
void MultiRadTrans::copySearchEntries(MultiRadTrans &tables)
{
    int j = 0;
    for(auto &i:tables.mOpacityTables)
    {
        mOpacityTables.at(j).searchQuantsOpacity = std::move(i.searchQuantsOpacity);
        j++;
    }
} 
// std::vector<LoadInTables>& MultiRadTrans::getSearchEntries()
// {
//     return mOpacityTables;
// }
void MultiRadTrans::calculateOpacities(GridData const *gridData, FixedData const &fixedData, bool skip_search)
{
    int oi = 0;
    double opamP, opamR;
    if(!skip_search)
    {
        mOpacityTables.at(0).findAllOPACITYIndicies(gridData);
        if(mMaterialInterface > 0)
        {
            mOpacityTables.at(1).findAllOPACITYIndicies(gridData);
        }
    }
    for(int i = 0; i < mNx; i++)
    {
        opamP = fixedData.opamP1;
        opamR = fixedData.opamR1;
        if((mMaterialInterface <= i) && (oi == 0))
        {
            oi = 1;
            opamP = fixedData.opamP2;
            opamR = fixedData.opamR2;
        }
        planckianOpacity(gridData, oi, i, opamP);
        rosslandOpacity(gridData, oi, i, opamR);
    }
}
void MultiRadTrans::planckianOpacity(GridData const *gridData, int const table, int const i, float const opamP)
{
    double eq11,eq12,eq21,eq22;
    int mNg = mOpacityTables.at(table).PhotonGrid.size() - 1;
    int elower_temp_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);
    int eupper_temp_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);

    int elower_density_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    int eupper_density_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    double opa_cap = 1e3;
    double min_cell_size = *std::min_element(gridData->Dx_k_1_2.begin(), gridData->Dx_k_1_2.end());
    double domain_size = gridData->CellWallCoord[mNx];

    //Format used by  tables 
    // int ncols = mOpacityTables.at(table).Density.size();
    // double eq11 = mOpacityTables.at(table).planckAbsorptionOpacity[mEnergyGroup + mNg*(elower_density_index + ncols * elower_temp_index)]; //q11
    // double eq12 = mOpacityTables.at(table).planckAbsorptionOpacity[mEnergyGroup + mNg*(elower_density_index + ncols * eupper_temp_index)]; //q12
    // double eq21 = mOpacityTables.at(table).planckAbsorptionOpacity[mEnergyGroup + mNg*(eupper_density_index + ncols * elower_temp_index)]; //q21
    // double eq22 = mOpacityTables.at(table).planckAbsorptionOpacity[mEnergyGroup + mNg*(eupper_density_index + ncols * eupper_temp_index)]; //q22
    //Format used by SPK Tables
    if(mNg == 1)
    {
        int ncols = mOpacityTables.at(table).OpacityTemperature.size();
        eq11 = mOpacityTables.at(table).planckAbsorptionOpacity[(elower_density_index * ncols + elower_temp_index)]; //q11
        eq12 = mOpacityTables.at(table).planckAbsorptionOpacity[(elower_density_index * ncols + eupper_temp_index)]; //q12
        eq21 = mOpacityTables.at(table).planckAbsorptionOpacity[(eupper_density_index * ncols + elower_temp_index)]; //q21
        eq22 = mOpacityTables.at(table).planckAbsorptionOpacity[(eupper_density_index * ncols + eupper_temp_index)]; //q22
    }
    else
    {
        //original shape: nNg, nRho,Nte = np.shape(Opacities)
        // Flat index formula = te + nTe*rho + mEnergy*nTe*nRho
        // Last part enters the right energy group table 
        // second part finds the right density row 
        // first part right Te 
        int nTe = mOpacityTables.at(table).OpacityTemperature.size();
        int nRho =mOpacityTables.at(table).OpacityDensity.size(); 
        eq11 = mOpacityTables.at(table).planckAbsorptionOpacity[elower_temp_index + nTe * elower_density_index + mEnergyGroup * nTe*nRho]; //q11
        eq12 = mOpacityTables.at(table).planckAbsorptionOpacity[eupper_temp_index + nTe * elower_density_index + mEnergyGroup * nTe*nRho]; //q12
        eq21 = mOpacityTables.at(table).planckAbsorptionOpacity[elower_temp_index + nTe * eupper_density_index + mEnergyGroup * nTe*nRho]; //q21
        eq22 = mOpacityTables.at(table).planckAbsorptionOpacity[eupper_temp_index + nTe * eupper_density_index + mEnergyGroup * nTe*nRho]; //q22
    }

    if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values) > gridData->TemperatureE[i])
    {
        mLocal.RadPlanckAbsorptionOpacity[i] = 0;
    }
    else if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values) > gridData->Density[i])
    {
        mLocal.RadPlanckAbsorptionOpacity[i] = 0;
    }
    else
    {
        mLocal.RadPlanckAbsorptionOpacity[i] = opamP * mOpacityTables.at(table).biLinearInterpolation(std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);
        if(mLocal.RadRossAbsorptionOpacity[i] < mLocalRossOpacityCap)
        {
            mLocal.RadRossAbsorptionOpacity[i] = mLocalRossOpacityCap;
        }
        //SpK unit conversion m^-1 -> m^2kg^-1

        //Gorgon Method
        // mLocal.RadPlanckAbsorptionOpacity[i] =  opamP * std::max(1./(domain_size*opa_cap),std::min(mLocal.RadPlanckAbsorptionOpacity[i],opa_cap/min_cell_size));
        mLocal.RadPlanckAbsorptionOpacity[i] /= gridData->Density[i];
        //Cap based methods
    }
    mLocal.RadPlanckEmissionOpacity[i] = mLocal.RadPlanckAbsorptionOpacity[i];
}
void MultiRadTrans::rosslandOpacity(GridData const *gridData,int const table, int const i, float const opamR)
{
    double eq11,eq12,eq21,eq22;
    int mNg = mOpacityTables.at(table).PhotonGrid.size() - 1;
    int elower_temp_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);
    int eupper_temp_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_index);

    int elower_density_index = std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    int eupper_density_index = std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_index);
    double opa_cap = 1e3;
    double min_cell_size = *std::min_element(gridData->Dx_k_1_2.begin(), gridData->Dx_k_1_2.end());
    double domain_size = gridData->CellWallCoord[mNx];

    ////Format used by  tables 
    // int ncols = mOpacityTables.at(table).Density.size();
    // double eq11 = mOpacityTables.at(table).rosslandOpacity[mEnergyGroup + mNg*(elower_density_index + ncols * elower_temp_index)]; //q11
    // double eq12 = mOpacityTables.at(table).rosslandOpacity[mEnergyGroup + mNg*(elower_density_index + ncols * eupper_temp_index)]; //q12
    // double eq21 = mOpacityTables.at(table).rosslandOpacity[mEnergyGroup + mNg*(eupper_density_index + ncols * elower_temp_index)]; //q21
    // double eq22 = mOpacityTables.at(table).rosslandOpacity[mEnergyGroup + mNg*(eupper_density_index + ncols * eupper_temp_index)]; //q22
    //Format used by SPK Tables
    if(mNg == 1)
    {
        int ncols = mOpacityTables.at(table).OpacityTemperature.size();
        eq11 = mOpacityTables.at(table).rosslandOpacity[(elower_density_index * ncols + elower_temp_index)]; //q11
        eq12 = mOpacityTables.at(table).rosslandOpacity[(elower_density_index * ncols + eupper_temp_index)]; //q12
        eq21 = mOpacityTables.at(table).rosslandOpacity[(eupper_density_index * ncols + elower_temp_index)]; //q21
        eq22 = mOpacityTables.at(table).rosslandOpacity[(eupper_density_index * ncols + eupper_temp_index)]; //q22
    }
    else
    {
        //original shape: nNg, nRho,Nte = np.shape(Opacities)
        // Flat index formula = te + nTe*rho + mEnergy*nTe*nRho
        // Last part enters the right energy group table 
        // second part finds the right density row 
        // first part right Te 
        int nTe = mOpacityTables.at(table).OpacityTemperature.size();
        int nRho =mOpacityTables.at(table).OpacityDensity.size(); 
        eq11 = mOpacityTables.at(table).rosslandOpacity[elower_temp_index + nTe * elower_density_index + mEnergyGroup * nTe*nRho]; //q11
        eq12 = mOpacityTables.at(table).rosslandOpacity[eupper_temp_index + nTe * elower_density_index + mEnergyGroup * nTe*nRho]; //q12
        eq21 = mOpacityTables.at(table).rosslandOpacity[elower_temp_index + nTe * eupper_density_index + mEnergyGroup * nTe*nRho]; //q21
        eq22 = mOpacityTables.at(table).rosslandOpacity[eupper_temp_index + nTe * eupper_density_index + mEnergyGroup * nTe*nRho]; //q22
    }

    if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values) > gridData->TemperatureE[i])
    {
        mLocal.RadRossAbsorptionOpacity[i] = 0;
    }
    else if (std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values) > gridData->Density[i])
    {
        mLocal.RadRossAbsorptionOpacity[i] = 0;
    }
    else
    {   
        mLocal.RadRossAbsorptionOpacity[i] =  opamR * mOpacityTables.at(table).biLinearInterpolation(std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].te_values),
                                                    std::get<0>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    std::get<1>(mOpacityTables.at(table).searchQuantsOpacity[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);
        if(mLocal.RadRossAbsorptionOpacity[i] < mLocalRossOpacityCap)
        {
            mLocal.RadRossAbsorptionOpacity[i] = mLocalRossOpacityCap;
        }
        //SpK unit conversion m^-1 -> m^2kg^-1
        mLocal.RadRossAbsorptionOpacity[i] /= gridData->Density[i];
        //GORGON method
        // mLocal.RadRossAbsorptionOpacity[i] = opamR * std::max(1./(domain_size*opa_cap),std::min(mLocal.RadRossAbsorptionOpacity[i],opa_cap/min_cell_size));
        //Floor Based Rescaling of opacity
    }
}
void MultiRadTrans::ficksDiffusionCoefficient(GridData const *gridData, FixedData const &fixedData)
{
    double larsen_limiter = fixedData.larsen_limiter;
    double diffusion, averaged_rossland_mean, rad_gradient, cell_wall_rad_energy, flux_limited_free_streaming;
    for(int i = 1; i < mNx; i++)
    {   //Convert from m^2kg^-1  -> m^-1
        averaged_rossland_mean = 0.5*(gridData->Density[i - 1]* mLocal.RadRossAbsorptionOpacity[i - 1] + 
                                    gridData->Density[i]* mLocal.RadRossAbsorptionOpacity[i]);
        //NOTE SQUARE ROOT FLUX LIMITER Morel, J. E. 2000, J. Quant. Spectrosc. Radiat. Transfer, 65, 769. 
        rad_gradient = (mLocal.RadEnergyDensity[i] - mLocal.RadEnergyDensity[i - 1]) / gridData->Dx_k[i - 1]; 
        cell_wall_rad_energy = (gridData->Dx_k_1_2[i - 1] * mLocal.RadEnergyDensity[i] +
                                       gridData->Dx_k_1_2[i] * mLocal.RadEnergyDensity[i - 1]) / 
                                       (gridData->Dx_k_1_2[i - 1] + gridData->Dx_k_1_2[i]);
                    

        // flux_limited_free_streaming = fixedData.radflxlm  * (pow(abs(rad_gradient), larsen_limiter) / pow(cell_wall_rad_energy, larsen_limiter));
        flux_limited_free_streaming = fixedData.radflxlm  * (abs(rad_gradient) * abs(rad_gradient)) / (cell_wall_rad_energy * cell_wall_rad_energy);

        if(averaged_rossland_mean < 1e-50)
        {
            diffusion =  0; //1e-10; // This will prevent Nans in the matrix. As long as this value is << dx_k_1_2 we are safe, if you go lower density make sure to reflect that here. 
        }
        else
        {
            double larsen = 1.0 / larsen_limiter;
            double opacity_larsen = pow(3.0 * averaged_rossland_mean, larsen_limiter);
            double bottom = pow(( opacity_larsen + flux_limited_free_streaming), larsen);

            diffusion = ((fixedData.SPEED_OF_LIGHT)) / bottom; 
        }
        mLocal.RadDiffusionCoefficient[i] = diffusion; 
    }
}
void MultiRadTrans::mLeftReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d)
{
    int i = 0;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    // mLeftBoundaryCol[0] = i ; mLeftBoundaryCol[1] = i + 1;
    std::get<0>(mLeftBoundaryValue) = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]));
    std::get<1>(mLeftBoundaryValue) = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
    mEnergyDensities[i] = mLocal.RadEnergyDensity[i];
    b[i] = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]));
    c[i] = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
    d[i] = mLocal.RadEnergyDensity[i]; 
}
void MultiRadTrans::mRightReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d)
{
    int i = mNx - 1;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    std::get<0>(mRightBoundaryValue) = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
    std::get<1>(mRightBoundaryValue) = 1 + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]));
    mEnergyDensities[i] = mLocal.RadEnergyDensity[i];
    a[i] = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1])); 
    b[i] = 1 + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1])); 
    d[i] = mLocal.RadEnergyDensity[i];
}
// void MultiRadTrans::mLeftVacuumSourceBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in, double *a, double *b, double *c, double*d)
// {
//     double D_0;
//     if(mLocal.RadRossAbsorptionOpacity[0] < 1e-50)
//     {
//         D_0 = 1e30;
//     }
//     else
//     {
//         D_0 = fixedData.SPEED_OF_LIGHT / (3 * gridData->Density[0]* mLocal.RadRossAbsorptionOpacity[0]); 
//     }
    
//     double D_eff =  1 / (1 + (4 / (fixedData.SPEED_OF_LIGHT * gridData->Dx_k_1_2[0])) * D_0);
//     double diff_1 = (timeData.Dt05 * mLocal.RadDiffusionCoefficient[1]) / (gridData->Dx_k_1_2[0] * gridData->Dx_k[0]); 
//     double diff_0 = (2*timeData.Dt05 * D_0 * D_eff) / (gridData->Dx_k_1_2[0] * gridData->Dx_k_1_2[0]); 
//     int i = 0;
//     double r = (timeData.Dt05 / gridData->Dx_k[i]);
//     std::get<0>(mLeftBoundaryValue) = 1 + diff_0 + diff_1 + timeData.Dt05; 
//     std::get<1>(mLeftBoundaryValue) =  -1 * diff_1;
//     mEnergyDensities[i] = mLocal.RadEnergyDensity[i] +  ((4 * diff_0) / (fixedData.SPEED_OF_LIGHT)) * f_in;
//     b[i] = 1 + diff_0 + diff_1 + timeData.Dt05; 
//     c[i] = -1 * diff_1;
//     d[i] = mLocal.RadEnergyDensity[i] +  ((4 * diff_0) / (fixedData.SPEED_OF_LIGHT)) * f_in; 
// }


void MultiRadTrans::ficksDiffusion(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData)
{   //Method from Bowers And Wilson p427 .. modified for my needs 
    //See lab book for currently implemented method. 
    //Matrix Equations being solved is:
    // (1 + r(Dhat_{k+1} + Dhat_{k}))U^{n+1}_{k+1/2} - rDhat_{k+1}U^{n+1}_{k+3/2} - rDhat_{k}U^{n+1}_{k-1/2} = U^{''}_{k+1/2}
    // r = dt_{n+1/2} / dr{k+1/2}
    // Dhat_{k+1} = D_r / dr{k+1}
    // D_r_k = (dr * SPEED_OF_LIGHT/3) / (rho_k * Kappa_Ross * dr_k + C_L mag(u_{k+1/2} - u_{k-1/2}) /0.5(u_{k+1/2} - u_{k-1/2}))
    // MultiRadTrans::createPETScObjects();
    std::vector<std::tuple<double, double, double>> diff_coef;
    bool direchletRight = false, direchletLeft = false;
    double val_0,val_1,val_2;
    diff_coef.resize(mNx);
    double a[mNx];
    double b[mNx];
    double c[mNx]; 
    double d[mNx];
    for(int i = 1; i < mNx - 1; i++)
    {
        double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
        val_0 = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
        val_1 = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]))
                    + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]));
        val_2 = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
        diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
        mEnergyDensities[i] = mLocal.RadEnergyDensity[i];
        a[i] = val_0;
        b[i] = val_1;
        c[i] = val_2;
        d[i] = mLocal.RadEnergyDensity[i];
    }
    if(fixedData.RightRadBoundaryCondition == "r")
    {
        mRightReflectiveBoundaryCondition(gridData, fixedData, timeData, a, b, c, d);
    }
    else
    {
        direchletRight = true;
        mEnergyDensities[mNx - 1] = fixedData.RightDirichletEnergyDensity;
        d[mNx - 1] =fixedData.RightDirichletEnergyDensity; 
        b[mNx - 1] = 1;
        a[mNx - 1] = 0;
    }
    if(fixedData.LeftRadBoundaryCondition == "r")
    {
        mLeftReflectiveBoundaryCondition(gridData, fixedData, timeData, a,b,c,d);
    }
    // else if((fixedData.LeftRadBoundaryCondition == "v") || (fixedData.LeftRadBoundaryCondition == "m"))
    // {
    //     double f_in = 0;
    //     if(fixedData.LeftRadBoundaryCondition  == "m")
    //     {
    //         f_in = 0.25*fixedData.SPEED_OF_LIGHT*fixedData.radiationConstant * pow(fixedData.LeftSourceTemperature, 4.0);
        // }
    //     mLeftVacuumSourceBoundaryCondition(gridData, fixedData, timeData, f_in,a,b,c,d);
    // }
    else
    {
        direchletLeft = true;
        mEnergyDensities[0] = fixedData.LeftDirichletEnergyDensity;
        d[0] = fixedData.LeftDirichletEnergyDensity; 
        b[0] = 1;
        c[0] = 0;
    }

    // matrixSolver.fillTriMatrix(diff_coef, mEnergyDensities, mLeftBoundaryValue, mRightBoundaryValue,
    //                              direchletLeft,direchletRight);
    // matrixSolver.solveMatrix();
    // matrixSolver.returnSolution(mLocal.RadEnergyDensity);
    mSolveTridiagonal(a, b, c, d, mNx);
    for(int i = 0; i < mNx; i ++)
    {
        mLocal.RadEnergyDensity[i] = d[i];
    }
    // std::cout << "Both calculated" << std::endl;
}

void MultiRadTrans::mImplicitRightReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d)
{
    int i = mNx - 1;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    std::get<0>(mRightBoundaryValue) = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
    std::get<1>(mRightBoundaryValue) = 1 + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]));
    mEnergyDensities[i] = mLocal.RadEnergyDensity[i];
    a[i] = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1])); 
    b[i] = 1 + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1])) 
             + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[i]* timeData.Dt05 * gridData->Density[i]; 
    d[i] = mLocal.RadEnergyDensity[i] + mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i];
}
void MultiRadTrans::mImplicitLeftReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d, double heating)
{
    int i = 0;
    double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
    // mLeftBoundaryCol[0] = i ; mLeftBoundaryCol[1] = i + 1;
    std::get<0>(mLeftBoundaryValue) = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]));
    std::get<1>(mLeftBoundaryValue) = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
    mEnergyDensities[i] = mLocal.RadEnergyDensity[i];
    b[i] = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]))
             + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[i]* timeData.Dt05 * gridData->Density[i];
    c[i] = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
    d[i] = mLocal.RadEnergyDensity[i]+ mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i] + timeData.Dt05 * heating;
    
}
void MultiRadTrans::mImplicitLeftVacuumSourceBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in, double *a, double *b, double *c, double*d)
{
    double D_0;
    if(mLocal.RadRossAbsorptionOpacity[0] < 1e-50)
    {
        D_0 = 1e30;
    }
    else
    {
        D_0 = fixedData.SPEED_OF_LIGHT / (3 * gridData->Density[0]* mLocal.RadRossAbsorptionOpacity[0]); 
    }
    
    int i = 0;
    double D_eff =  1 / (1 + (4 / (fixedData.SPEED_OF_LIGHT * gridData->Dx_k_1_2[0])) * D_0);
    double diff_1 = (timeData.Dt05 * mLocal.RadDiffusionCoefficient[1]) / (gridData->Dx_k_1_2[0] * gridData->Dx_k[0]); 
    // double diff_0 = (2*timeData.Dt05 * D_0 * D_eff) / (gridData->Dx_k_1_2[0]); 
    double diff_0 = (D_0*D_eff * 2 * timeData.Dt05) / (gridData->Dx_k_1_2[0] * gridData->Dx_k_1_2[0]);
    b[i] = 1 + diff_0 + diff_1 + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[0]* timeData.Dt05 * gridData->Density[i]; 
    c[i] = -1 * diff_1;
    // d[i] = mLocal.RadEnergyDensity[i] +  ((4 * diff_0) / (fixedData.SPEED_OF_LIGHT)) * f_in + mLocal.RadFFEmission[0]* timeData.Dt05 * gridData->Density[i]; 
    // d[i] = mLocal.RadEnergyDensity[0] + (timeData.Dt05 * 8 * f_in * D_eff * D_0) / (fixedData.SPEED_OF_LIGHT*gridData->Dx_k_1_2[0] * gridData->Dx_k_1_2[0]) + mLocal.RadFFEmission[0]* timeData.Dt05 * gridData->Density[0];
    d[i] = mLocal.RadEnergyDensity[0] + mLocal.RadFFEmission[0]* timeData.Dt05 * gridData->Density[0];
}
void MultiRadTrans::fullyImplicit(GridData const* gridData, FixedData const& fixedData, TimeStep const &timeData, std::vector<double> heating)
{
    //See lab book for currently implemented method. 
    //Matrix Equations being solved is:
    // (1 + r(Dhat_{k+1} + Dhat_{k}))U^{n+1}_{k+1/2} - rDhat_{k+1}U^{n+1}_{k+3/2} - rDhat_{k}U^{n+1}_{k-1/2} = U^{''}_{k+1/2}
    // r = dt_{n+1/2} / dr{k+1/2}
    // Dhat_{k+1} = D_r / dr{k+1}
    // D_r_k = (dr * SPEED_OF_LIGHT/3) / (rho_k * Kappa_Ross * dr_k + C_L mag(u_{k+1/2} - u_{k-1/2}) /0.5(u_{k+1/2} - u_{k-1/2}))
    // MultiRadTrans::createPETScObjects();
    // ficksDiffusionCoefficient(gridData, fixedData);
    // calcFreeFreeEmission(gridData, fixedData);
    std::vector<std::tuple<double, double, double>> diff_coef;
    bool direchletRight = false, direchletLeft = false;
    double val_0,val_1,val_2;
    diff_coef.resize(mNx);
    double a[mNx];
    double b[mNx];
    double c[mNx]; 
    double d[mNx];
    for(int i = 1; i < mNx - 1; i++)
    {
        double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
        val_0 = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
        val_1 = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]))
                    + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]))
                    + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[i]* timeData.Dt05 * gridData->Density[i];
        val_2 = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
        // diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
        // mEnergyDensities[i] = mLocal.RadEnergyDensity[i] + mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i];
        a[i] = val_0;
        b[i] = val_1;
        c[i] = val_2;
        d[i] = mLocal.RadEnergyDensity[i] + mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i] + timeData.Dt05 * heating[i];
    }
    if(fixedData.RightRadBoundaryCondition == "r")
    {
        mImplicitRightReflectiveBoundaryCondition(gridData, fixedData, timeData, a, b, c, d);
    }
    else
    {
        direchletRight = true;
        mEnergyDensities[mNx - 1] = fixedData.RightDirichletEnergyDensity;
        d[mNx - 1] =fixedData.RightDirichletEnergyDensity; 
        b[mNx - 1] = 1;
        a[mNx - 1] = 0;
    }
    if(fixedData.LeftRadBoundaryCondition == "r")
    {
        mImplicitLeftReflectiveBoundaryCondition(gridData, fixedData, timeData, a,b,c,d, heating[0]);
    }
    else if((fixedData.LeftRadBoundaryCondition == "v") || (fixedData.LeftRadBoundaryCondition == "m"))
    {
        double f_in = 0;
        if(fixedData.LeftRadBoundaryCondition  == "m")
        {
            f_in = fixedData.radiationConstant * pow(fixedData.LeftSourceTemperature, 4.0);
        }

        mImplicitLeftVacuumSourceBoundaryCondition(gridData, fixedData, timeData, f_in,a,b,c,d);
    }
    else
    {
        direchletLeft = true;
        mEnergyDensities[0] = fixedData.LeftDirichletEnergyDensity;
        d[0] = fixedData.LeftDirichletEnergyDensity; 
        b[0] = 1;
        c[0] = 0;
    }

    mSolveTridiagonal(a, b, c, d, mNx);
    for(int i = 0; i < mNx; i ++)
    {
        mLocal.RadEnergyDensity[i] = d[i];
    }

}
void MultiRadTrans::fullyImplicit(GridData const* gridData, FixedData const& fixedData, TimeStep const &timeData)
{
    //See lab book for currently implemented method. 
    //Matrix Equations being solved is:
    // (1 + r(Dhat_{k+1} + Dhat_{k}))U^{n+1}_{k+1/2} - rDhat_{k+1}U^{n+1}_{k+3/2} - rDhat_{k}U^{n+1}_{k-1/2} = U^{''}_{k+1/2}
    // r = dt_{n+1/2} / dr{k+1/2}
    // Dhat_{k+1} = D_r / dr{k+1}
    // D_r_k = (dr * SPEED_OF_LIGHT/3) / (rho_k * Kappa_Ross * dr_k + C_L mag(u_{k+1/2} - u_{k-1/2}) /0.5(u_{k+1/2} - u_{k-1/2}))
    // MultiRadTrans::createPETScObjects();
    // ficksDiffusionCoefficient(gridData, fixedData);
    // calcFreeFreeEmission(gridData, fixedData);
    if(mNx > 1) //REF how would diffusion even work in a single cell - just boundary conditions? 
    {
        std::vector<std::tuple<double, double, double>> diff_coef;
        bool direchletRight = false, direchletLeft = false;
        double val_0,val_1,val_2;
        diff_coef.resize(mNx);
        double a[mNx];
        double b[mNx];
        double c[mNx]; 
        double d[mNx];
        for(int i = 1; i < mNx - 1; i++)
        {
            double r = (timeData.Dt05 / gridData->Dx_k_1_2[i]);
            val_0 = -1*r * (mLocal.RadDiffusionCoefficient[i]) / ((gridData->Dx_k[i - 1]));
            val_1 = 1 + r * mLocal.RadDiffusionCoefficient[i + 1]/((gridData->Dx_k[i]))
                        + r * mLocal.RadDiffusionCoefficient[i] / ((gridData->Dx_k[i - 1]))
                        + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[i]* timeData.Dt05 * gridData->Density[i]
                        + mLeakSourceCoefficient[i];
            val_2 = -1*r * (mLocal.RadDiffusionCoefficient[i + 1]) /((gridData->Dx_k[i]));
            // diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
            // mEnergyDensities[i] = mLocal.RadEnergyDensity[i] + mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i];
            a[i] = val_0;
            b[i] = val_1;
            c[i] = val_2;
            d[i] = mLocal.RadEnergyDensity[i] + mLocal.RadFFEmission[i] * timeData.Dt05 * gridData->Density[i];
            
            
        }
        if(fixedData.RightRadBoundaryCondition == "r")
        {
            mImplicitRightReflectiveBoundaryCondition(gridData, fixedData, timeData, a, b, c, d);
        }
        else if(fixedData.RightRadBoundaryCondition == "c")
        {
            b[mNx - 1] = 1 + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[mNx - 1]* timeData.Dt05 * gridData->Density[mNx - 1];
            d[mNx - 1] = mLocal.RadEnergyDensity[mNx - 1] + mLocal.RadFFEmission[mNx - 1] * timeData.Dt05 * gridData->Density[mNx - 1];
            a[mNx -1 ] = 0;
        }
        else
        {
            direchletRight = true;
            mEnergyDensities[mNx - 1] = fixedData.RightDirichletEnergyDensity;
            d[mNx - 1] =fixedData.RightDirichletEnergyDensity; 
            b[mNx - 1] = 1;
            a[mNx - 1] = 0;
        }
        if(fixedData.LeftRadBoundaryCondition == "r")
        {
            mImplicitLeftReflectiveBoundaryCondition(gridData, fixedData, timeData, a,b,c,d, 0.0);
        }
        else if (fixedData.LeftRadBoundaryCondition == "c")
        {
            b[0] = 1 + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[0]* timeData.Dt05 * gridData->Density[0];;
            c[0] = 0;
            d[0] = mLocal.RadEnergyDensity[0] + mLocal.RadFFEmission[0] * timeData.Dt05 * gridData->Density[0];

        }
        else if((fixedData.LeftRadBoundaryCondition == "v") || (fixedData.LeftRadBoundaryCondition == "m"))
        {
            double f_in = 0;
            if(fixedData.LeftRadBoundaryCondition  == "m")
            {
                f_in = fixedData.radiationConstant * pow(fixedData.LeftSourceTemperature, 4.0);
            }

            mImplicitLeftVacuumSourceBoundaryCondition(gridData, fixedData, timeData, f_in,a,b,c,d);
        }
        else
        {
            direchletLeft = true;
            mEnergyDensities[0] = fixedData.LeftDirichletEnergyDensity;
            d[0] = fixedData.LeftDirichletEnergyDensity; 
            b[0] = 1;
            c[0] = 0;
        }

        mSolveTridiagonal(a, b, c, d, mNx);
        for(int i = 0; i < mNx; i ++)
        {
            if(d[i] < mLocalEnergyDensityCap) //REF We hard capping the energy density
            {
                assert(d[i] >= 0);
                mLocal.RadEnergyDensity[i] = mLocalEnergyDensityCap;
            }
            else
            {
                mLocal.RadEnergyDensity[i] = d[i];
            }
        }
    }
    else
    {
        mLocal.RadEnergyDensity[0] = (timeData.Dt05 * gridData->Density[0] * mLocal.RadFFEmission[0] + mLocal.RadEnergyDensity[0]) / 
                                        (1 + fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[0] * timeData.Dt05 * gridData->Density[0]);
    }
}

void MultiRadTrans::calcFreeFreeAbsorption(GridData const *gridDataT,FixedData const&fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        mLocal.RadFFAbsorb[i] = fixedData.SPEED_OF_LIGHT * mLocal.RadPlanckAbsorptionOpacity[i] * mLocal.RadEnergyDensity[i]; 
        assert(mLocal.RadFFAbsorb[i]>=0);
    }
}
double MultiRadTrans::mPlanckPiLFunction(double x, int l, bool without_1)
{
    //cite Bradley Clark J.Comp Phys 70 1987
    auto series = [](double x, int l, int power)
    {double sum = 0;
    for(int i = 1; i <= l; i++)
    {
        sum += exp(-1*i*x)/pow(i,power); // Orig Implementation
    }
    return sum;};

    double pi_l;
    if(!without_1)
    {
        // pi_l = 1 + (15.0/pow(M_PI,4))* (-1.0 * pow(x,3)*series(x,l,1) - 3.0 * pow(x,2)*series(x,l,2) - 6.0 * x *series(x,l,3) - 6.0 * series(x,l,4));
        // pi_l = 1 + (15.0/pow(M_PI,4))* (-1.0 * pow(x,3)*series(x,l,1) - 3.0 * pow(x,2)*series(x,l,2) - 6.0 * x *series(x,l,3) - 6.0 * series(x,l,4));
        pi_l = 1 + (15.0/(M_PI *M_PI * M_PI * M_PI)) * (-1.0 * x*x*x * series(x,l,1) - 3.0*x*x*series(x,l,2) - 6.0 *x *series(x,l,3) - 6.0 * series(x,l,4));
    }
    else
    {
        // pi_l = (15.0/pow(M_PI,4))* (-1.0 * pow(x,3)*series(x,l,1) - 3.0 * pow(x,2)*series(x,l,2) - 6.0 * x *series(x,l,3) - 6.0 * series(x,l,4));
        // pi_l = (15.0/pow(M_PI,4))* (-1.0 * pow(x,3)*series(x,l,1) - 3.0 * pow(x,2)*series(x,l,2) - 6.0 * x *series(x,l,3) - 6.0 * series(x,l,4));
        pi_l = (15.0/(M_PI *M_PI * M_PI * M_PI)) * (-1.0 * x*x*x * series(x,l,1) - 3.0*x*x*series(x,l,2) - 6.0 *x *series(x,l,3) - 6.0 * series(x,l,4));
    }
    
    return pi_l;
}
double MultiRadTrans::mReturnGamma(double x)
{
    //cite Bradley Clark J.Comp Phys 70 1987
    // auto pi_n = [](auto x){return ((15.0/pow(M_PI,4))*((1.0/3.0)*pow(x,3) - (1.0/8.0)*pow(x,4) +(1.0/60.0)*pow(x,5) - (1.0/5040.0) * pow(x,7) +(1.0/272160.0)*pow(x,9)));};//hardcoded cutoff at 9 
    auto pi_n = [](auto x){return ((15.0/(M_PI *M_PI * M_PI * M_PI))*((1.0/3.0)*x*x*x - (1.0/8.0)*x*x*x*x +(1.0/60.0)*x*x*x*x*x - (1.0/5040.0) * x*x*x*x*x*x*x +(1.0/272160.0)*x*x*x*x*x*x*x*x*x));};//hardcoded cutoff at 9 
    if(x <= 1e-30)
    {
        return 0.0;
    }
    if(x > 1e30)
    {
        return 1.0;
    }
    
    double pi_n_hat = pi_n(x);
    double pi_l = mPlanckPiLFunction(x, 3, false); //REF hardcoding 3 here recommended config in paper is 9,3 for efficiency vs accuracy tradeoff
    if (pi_n_hat < pi_l)
    {
        return pi_n_hat;
    }
    else
    {
        double pi_l_1 = mPlanckPiLFunction(x, 3, true); //REF hardcoding 3 here recommended config in paper is 9,3 for efficiency vs accuracy tradeoff
        return pi_l_1;
    }
}
void MultiRadTrans::planckIntegral(GridData const *gridDataT, FixedData const&fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        double x1 = mOpacityTables.at(0).PhotonGrid[mEnergyGroup] / (gridDataT->TemperatureE[i] * fixedData.k_to_ev);
        double x2 = mOpacityTables.at(0).PhotonGrid[mEnergyGroup + 1] / (gridDataT->TemperatureE[i] * fixedData.k_to_ev);
        double tmp0, tmp1;
        // tmp.push_back(mReturnGamma(x));    
        tmp0 = mReturnGamma(x1);
        tmp1 = mReturnGamma(x2);
        double diff = tmp1 - tmp0;
        if(diff < 0)
        {
            diff = diff + 1;
        }
        // mPlanck[i] = (fixedData.STEFAN_BOLTZMANN_CONSTANT/(M_PI)) * pow(gridDataT->TemperatureE[i], 4.0) * diff; 
        mPlanck[i] = (fixedData.STEFAN_BOLTZMANN_CONSTANT/(M_PI)) * gridDataT->TemperatureE[i]*gridDataT->TemperatureE[i]*gridDataT->TemperatureE[i]*gridDataT->TemperatureE[i] * diff; 
        mPlanckIntegral[i] = diff;
    }
}

void MultiRadTrans::calcFreeFreeEmission(GridData const *gridDataT,FixedData const&fixedData)
{
    //double constant = 3.492717789885802e-08; //8*pi*k^4/c^2*h^3
    for(int i = 0; i < mNx; i++)
    {
        double planckian = mPlanck[i];
        mLocal.RadFFEmission[i] = 4.0 * M_PI * mLocal.RadPlanckEmissionOpacity[i]* planckian;  
        assert(mLocal.RadFFEmission[i]>=0);
        if (mLocal.RadFFEmission[i] < mLocalEmissionCap) //REF HARD CAP VALue
        {
            mLocal.RadFFEmission[i]  = mLocalEmissionCap;
        }
    }
}
void MultiRadTrans::updateRadiationEnergyDensity(GridData const *gridDataT,FixedData const&fixedData, TimeStep const &timeData)
{      
    // Update for emission
    // Emission(J) - Absorption(A) WKg-1
    // U'' = U' + rho*(4 *sigma* kpE(Te) * Te^4 - ckpA(Te, Tr) *Ur)
    calcFreeFreeAbsorption(gridDataT, fixedData);
    calcFreeFreeEmission(gridDataT, fixedData);

    for(int i = 0; i < mNx; i++)
    {
        mLocal.RadEnergyDensity[i] += gridDataT->Density[i] * timeData.Dt05 * (mLocal.RadFFEmission[i] - mLocal.RadFFAbsorb[i]);
    }
}

void MultiRadTrans::updateRadiationQuantities(GridData const *gridDataT, FixedData const &fixedData, TimeStep const &timeData, int updateType)
{
    if(updateType == 1)
    {   /*Update for emission
        Update Chain
        Calculate Rad Coefficient
        Calculate new heat Capacity which accounts for the radiation emission. 
        Calculate Tr'' 
        Update U''
        U'' = U' + J - A*/
        calculateOpacities(gridDataT, fixedData, true);
        planckIntegral(gridDataT, fixedData);
        updateRadiationEnergyDensity(gridDataT, fixedData, timeData);
        // LeakSource(fixedData, timeData);
        // updateRadiationTemperature(gridDataT, fixedData);
    }
    else if(updateType == 2)
    {   //Update for Ficks Diffusion
        //Update Chain
        //Calculate opacities
        //Calcualte Diffusion coefficient
        //Calculate the diffusion to find U^n+1
        //U^(n+1) = U_ficks
        calculateOpacities(gridDataT, fixedData, false);
        // RosslandFreeFree(gridDataT, fixedData);
        double sum_check = std::accumulate(mLocal.RadRossAbsorptionOpacity.begin(), mLocal.RadRossAbsorptionOpacity.end(), 0.0);
        if(sum_check > 1e-50)
        {
            ficksDiffusionCoefficient(gridDataT, fixedData);
            ficksDiffusion(gridDataT, fixedData, timeData);
        }
    }
}
void MultiRadTrans::updateRadiationQuantities(GridData const *gridDataT, FixedData const &fixedData, TimeStep const &timeData)
{
        if(mEnergyGroup == 0)
        {
            calculateOpacities(gridDataT, fixedData, false);
        }
        else
        {
            calculateOpacities(gridDataT, fixedData, true);
        }
        // setAnalyticalOpacity(gridDataT, fixedData);
        planckIntegral(gridDataT, fixedData);
        calcFreeFreeEmission(gridDataT, fixedData);
        double sum_check = std::accumulate(mLocal.RadRossAbsorptionOpacity.begin(), mLocal.RadRossAbsorptionOpacity.end(), 0.0);
        if((fixedData.leakArea1 != 0)||(fixedData.leakArea2 != 0))
        {
            LeakSource(gridDataT, fixedData, timeData);
        }
        if(sum_check > 1e-50)
        {
            std::vector<double> heating(fixedData.Nx);
            ficksDiffusionCoefficient(gridDataT, fixedData);
            fullyImplicit(gridDataT, fixedData, timeData);
        }
        calcFreeFreeAbsorption(gridDataT, fixedData);
}
void MultiRadTrans::updateRadiationTemperature(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadTemperature[i] = pow(((fixedData.SPEED_OF_LIGHT /  (4 * fixedData.STEFAN_BOLTZMANN_CONSTANT)) * gridData->RadEnergyDensity[i]), 0.25);
    }
}
void MultiRadTrans::setAnalyticalOpacity(GridData const *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        if(gridData->TemperatureE[i] > 1e3*fixedData.ev_to_k)
        {
            mLocal.RadRossAbsorptionOpacity[i] =  1e-1 * 6.0e3* (fixedData.Z[i] / fixedData.Ar[i]) *(gridData->TemperatureE[i] * 1e-3*fixedData.k_to_ev);
        }
        else
        {
            if(fixedData.Ar[i] == 64)
            {
                mLocal.RadRossAbsorptionOpacity[i] = 1e-1*13.9 * pow(gridData->Density[i]*1e-6, 0.29) / pow(gridData->TemperatureE[i]/11594000, 2.21);
                mLocal.RadPlanckAbsorptionOpacity[i] =1e-1*100.81 * pow(gridData->Density[i]*1e-6, 0.355) / pow(gridData->TemperatureE[i]/11594000, 2.12);
            }
            else
            {
                mLocal.RadRossAbsorptionOpacity[i] = 1e-1*1.9e-1 * pow(gridData->Density[i]*1e-3, 0.69) / pow(gridData->TemperatureE[i]*1e-3*fixedData.k_to_ev, 2.81);
                mLocal.RadPlanckAbsorptionOpacity[i] =1e-1*1.81e-1 * pow(gridData->Density[i]*1e-3, 0.755) / pow(gridData->TemperatureE[i]*1e-3*fixedData.k_to_ev, 2.52);

            }
        }
        mLocal.RadPlanckEmissionOpacity[i] = mLocal.RadRossAbsorptionOpacity[i];
    }
}
void MultiRadTrans::gather(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadFFEmission[i] += mLocal.RadFFEmission[i] ;
        gridData->RadFFAbsorb[i] += mLocal.RadFFAbsorb[i];
        gridData->RadEnergyDensity[i] += mLocal.RadEnergyDensity[i];
        if(gridData->ALlRadEnergyDensity.size() == 0)
        {
            gridData->RadRossAbsorptionOpacity.resize(mNx*mNg);
            gridData->RadPlanckEmissionOpacity.resize(mNx*mNg);
            gridData->RadPlanckAbsorptionOpacity.resize(mNx*mNg);
            gridData->ALlRadEnergyDensity.resize(mNx*mNg);
        }
        if(gridData->ALlRadEnergyDensity.size() > 0)
        {
            gridData->ALlRadEnergyDensity[mEnergyGroup*mNx + i] = mLocal.RadEnergyDensity[i];
            gridData->RadRossAbsorptionOpacity[mEnergyGroup*mNx + i] = mLocal.RadRossAbsorptionOpacity[i];
            gridData->RadPlanckEmissionOpacity[mEnergyGroup*mNx + i] = mLocal.RadPlanckEmissionOpacity[i];
            gridData->RadPlanckAbsorptionOpacity[mEnergyGroup*mNx + i] = mLocal.RadPlanckAbsorptionOpacity[i];
        }
    }
}
void MultiRadTrans::setPlanckianOpacity(double value)
{
    std::fill(mLocal.RadPlanckEmissionOpacity.begin(), mLocal.RadPlanckEmissionOpacity.end(),value); 
    std::fill(mLocal.RadPlanckAbsorptionOpacity.begin(), mLocal.RadPlanckAbsorptionOpacity.end(),value); 
}
void MultiRadTrans::setRossOpacity( double value)
{
    std::fill(mLocal.RadRossAbsorptionOpacity.begin(), mLocal.RadRossAbsorptionOpacity.end(), value);
}
void MultiRadTrans::setRadEnergyDensity(double value)
{
    std::fill(mLocal.RadEnergyDensity.begin(), mLocal.RadEnergyDensity.end(), value); 
}
void MultiRadTrans::setIntegratedRadEnergyDensity(std::vector<double> totalRadEnergyDensity)
{
    for(int i = 0; i < mNx; i++)
    {
        mLocal.RadEnergyDensity[i] = ((mOpacityTables.at(0).PhotonGrid[mEnergyGroup] - mOpacityTables.at(0).PhotonGrid[mEnergyGroup - 1]) / 2) * totalRadEnergyDensity[i]; 
    }
}
void MultiRadTrans::setRadEnergyDensity( std::vector<double> values)
{
    for(int i = 0; i < values.size(); i++)
    {
        mLocal.RadEnergyDensity[i] = values[i];
    }
}
void MultiRadTrans::setTotalEnergyGroup(int totalGroups)
{
    mNg = totalGroups;
}
void MultiRadTrans::setPlanckIntegral(std::vector<double> &values)
{
    int i = 0;
    for(auto val: values)
    {
        mPlanck[i] = val;
        i++; 
    }
}
void MultiRadTrans::setDiffusionCoefficient(std::vector<double> &values)
{
    int i = 0;
    for(auto val: values)
    {
        mLocal.RadDiffusionCoefficient[i] = val;
        i++; 
    }
}
std::vector<double> MultiRadTrans::getLocalEnergyDensity()
{
    return mLocal.RadEnergyDensity;
}
void MultiRadTrans::heatRegion(std::vector<double> heating,TimeStep const &timeData)
{
    for(int i = 0; i < mNx; i++)
    {
        mLocal.RadEnergyDensity[i] += timeData.Dt05 * heating[i];
    }
}
std::vector<double> MultiRadTrans::getPlanckIntegral()
{
    return mPlanck;
}
// void MultiRadTrans::LeakSource(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData)
// {
//     double summed_volume = 0;
//     for(int i = 0; i < mNx; i++)
//     {
//         if(mLeakAreaVector[i] > 0)
        // {
        //     summed_volume += (1.0 / gridData->Dx_k_1_2[i]);
    //     }
    // }
    // for(int i = 0; i < mNx; i++ )
    // {
        // mLocal.RadEnergyDensity[i] -=  0.25 * mLeakAreaVector[i] * mLocal.RadEnergyDensity[i] * timeData.Dt05 * fixedData.SPEED_OF_LIGHT * summed_volume;
//     }
// }
void MultiRadTrans::LeakSource(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData)
{
    double summed_volume = 0;
    double summed_area = 0;
    for(int i = 0; i < mNx; i++)
    {
        if(mLeakAreaVector[i] > 0)
        {
            summed_area+= (gridData->Dx_k_1_2[i]);
        }
    }
    summed_volume = 1/(1e-4 * summed_area); //assume a 1cm 1cm box hence 1e-4 area term
    for(int i = 0; i < mNx; i++ )
    {
        mLeakSourceCoefficient[i] =  0.25 * mLeakAreaVector[i] * timeData.Dt05 * fixedData.SPEED_OF_LIGHT * summed_volume;
    }
}

void MultiRadTrans::mSolveTridiagonal(double* a, double* b, double* c, double* d, int n) {
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

void MultiRadTrans::getOpacities(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->RadRossAbsorptionOpacity[i] = mLocal.RadRossAbsorptionOpacity[i];
        gridData->RadPlanckEmissionOpacity[i] = mLocal.RadPlanckEmissionOpacity[i];
        gridData->RadPlanckAbsorptionOpacity[i] = mLocal.RadPlanckEmissionOpacity[i];
    }
}