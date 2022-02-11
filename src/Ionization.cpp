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
#include "HyKiCT/Ionization.h"

Ionization::Ionization(int startNx, int nx)
{
   mStartNx = startNx;
   mNx = nx;
}
Ionization::Ionization(int startNx, int nx, std::unique_ptr<LoadInTables> &table)
{
   mStartNx = startNx;
   mNx = nx;
   mIonizationTable = std::move(table);
}
void Ionization::mFEOSUpdateIonisation(GridData *gridData)
{
    for(int i = mStartNx; i < mNx; i++)
    {
        // std::vector<double> e_feos_values_and_bounds = mReturnFeosValues(mSearchIndiciesElectron[i], mIonisationFeosTable, i * 4); 
        //Recall everything is packed as 4s hence *4
        int lower_temp_index = std::get<0>(mIonizationTable->searchQuantsElectron[i].te_index);
        int upper_temp_index = std::get<1>(mIonizationTable->searchQuantsElectron[i].te_index);

        int lower_density_index = std::get<0>(mIonizationTable->searchQuantsElectron[i].rho_index);
        int upper_density_index = std::get<1>(mIonizationTable->searchQuantsElectron[i].rho_index);
        int ncols = mIonizationTable->FEOSTemperature.size();
        double eq11 = mIonizationTable->ionisationFeosTable[lower_density_index * ncols + lower_temp_index]; //q11
        double eq12 = mIonizationTable->ionisationFeosTable[lower_density_index * ncols + upper_temp_index]; //q12
        double eq21 = mIonizationTable->ionisationFeosTable[upper_density_index * ncols + lower_temp_index]; //q21
        double eq22  = mIonizationTable->ionisationFeosTable[upper_density_index * ncols + upper_temp_index]; //q22
        
        // double elower_temp = e_feos_values_and_bounds[0] * fixedData.ev_to_k;
        // double eupper_temp = e_feos_values_and_bounds[1] * fixedData.ev_to_k;
        // double elower_density = e_feos_values_and_bounds[2];
        // double eupper_density = e_feos_values_and_bounds[3];
        // double eq11 = e_feos_values_and_bounds[4];
        // double eq12 = e_feos_values_and_bounds[5];
        // double eq21 = e_feos_values_and_bounds[6];
        // double eq22 = e_feos_values_and_bounds[7];
        

        gridData->Zbar[i] = mIonizationTable->biLinearInterpolation(std::get<0>(mIonizationTable->searchQuantsElectron[i].te_values),
                                                    std::get<1>(mIonizationTable->searchQuantsElectron[i].te_values),
                                                    std::get<0>(mIonizationTable->searchQuantsElectron[i].rho_values),
                                                    std::get<1>(mIonizationTable->searchQuantsElectron[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);

        // gridData->Zbar[i] = mBiLinearInterpolation(elower_temp, eupper_temp, elower_density, eupper_density,
        //                          eq11, eq12, eq21, eq22, gridData->TemperatureE[i], gridData->Density[i]);
    }
}
void Ionization::updateIonization(GridData *gridData, FixedData const &fixedData,Switches const &switchData)
{
    if(switchData.FullyIonized)
    {
        Ionization::mFullyIonized(gridData, fixedData);
    }
    else if(switchData.IonizationTables)
    {
        mIonizationTable->findAllFeosIndicies(gridData);
        mFEOSUpdateIonisation(gridData);
    }
    else
    {
        Ionization::mAproxThomasFermiIonize(gridData, fixedData);
    }
}
void Ionization::mAproxThomasFermiIonize(GridData *gridData, FixedData const &fixedData)
{
    //Aprox thomas fermi ionisation
    //Constants from Rage simulations.
    float a1 = 0.003323467;
    float a2 = 0.97183224;
    float a3 = 9.26148E-5;
    float a4 = 3.1016524;
    float b0 = -1.762999;
    float b1 = 1.4317567;
    float b2 = 0.31546338;
    float c1 = -0.36666667;
    float c2 = 0.98333333;
    float alpha = 14.3139316;
    float beta = 0.66240046;
    float fthrds = 1.33333333;
    double convertKtoEv = fixedData.BOLTZMANN_CONSTANT / fixedData.ELECTRON_CHARGE;
    double convertMtoCm = 1e-6;

    for(int i = 0; i < mNx; i++)
    {
        double tzero = (gridData->TemperatureE[i] * convertKtoEv) / pow(fixedData.Z[i],fthrds);
        
        double r = (gridData->NumberDensityI[i] * convertMtoCm) / (fixedData.Z[i]);
        double tf = tzero / (1.0 + tzero);

        double aa = (a1* pow(tzero,a2)) + (a3 * pow(tzero,a4));
        double b = -exp(b0 + (b1 * tf) + (b2 * pow(tf,7)));
        double SPEED_OF_LIGHT = c2 + (c1*tf);

        double q1 = aa * pow(r, b);
        double q = pow(r, SPEED_OF_LIGHT) + pow(q1, SPEED_OF_LIGHT);
        double cm = 1 / SPEED_OF_LIGHT;
        double q_redef = pow(q, cm);

        double x = alpha * pow(q_redef, beta);

        double fraction = x / (1.0 + x + sqrt(1.0 + (2.0 * x)));
        gridData->Zbar[i] = fraction * fixedData.Z[i];
    }
}

void Ionization::mFullyIonized(GridData *gridData, FixedData const &fixedData)
{
    for(int i = mStartNx; i < mNx; i++)
    {
        gridData ->Zbar[i] = fixedData.Z[i];
    }
}