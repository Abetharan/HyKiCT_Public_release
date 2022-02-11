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
#include "HyKiCT/LoadInEoS.h"

LoadInEoS::LoadInEoS(int startNx, int nx, std::unique_ptr<LoadInTables> &tables)
{
    mStartNx = startNx;
    mNx = nx;
    mTables = std::move(tables); 
}

void LoadInEoS::FindIndicies(GridData *gridData)
{
    mTables->findAllFeosIndicies(gridData);
}
void LoadInEoS::setSearchIndicies(std::unique_ptr<LoadInTables> &table)
{
    mTables->searchQuantsElectron = std::move(table->searchQuantsElectron);
    mTables->searchQuantsIon = std::move(table->searchQuantsIon);
}
void LoadInEoS::updatePressureTerms(GridData *gridData)
{  

    for(int i = mStartNx; i < mNx; i++)
    {
        
        int ncols = mTables->FEOSTemperature.size();
        int elower_temp_index = std::get<0>(mTables->searchQuantsElectron[i].te_index);
        int eupper_temp_index = std::get<1>(mTables->searchQuantsElectron[i].te_index);

        int elower_density_index = std::get<0>(mTables->searchQuantsElectron[i].rho_index);
        int eupper_density_index = std::get<1>(mTables->searchQuantsElectron[i].rho_index);
        
        double eq11 = mPreMultiE*mTables->pressureElectronFeosTable[elower_density_index * ncols + elower_temp_index]; //q11
        double eq12 = mPreMultiE*mTables->pressureElectronFeosTable[elower_density_index * ncols + eupper_temp_index]; //q12
        double eq21 = mPreMultiE*mTables->pressureElectronFeosTable[eupper_density_index * ncols + elower_temp_index]; //q21
        double eq22 = mPreMultiE*mTables->pressureElectronFeosTable[eupper_density_index * ncols + eupper_temp_index]; //q22

        int ilower_temp_index = std::get<0>(mTables->searchQuantsIon[i].te_index);
        int iupper_temp_index = std::get<1>(mTables->searchQuantsIon[i].te_index);

        int ilower_density_index = std::get<0>(mTables->searchQuantsIon[i].rho_index);
        int iupper_density_index = std::get<1>(mTables->searchQuantsIon[i].rho_index);

        double iq11 = mPreMultiI*mTables->pressureIonFeosTable[ilower_density_index * ncols + ilower_temp_index]; //q11
        double iq12 = mPreMultiI*mTables->pressureIonFeosTable[ilower_density_index * ncols + iupper_temp_index]; //q12
        double iq21 = mPreMultiI*mTables->pressureIonFeosTable[iupper_density_index * ncols + ilower_temp_index]; //q21
        double iq22  = mPreMultiI*mTables->pressureIonFeosTable[iupper_density_index * ncols + iupper_temp_index]; //q22
        
        gridData->PressureE[i] = mTables->biLinearInterpolation(std::get<0>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<0>(mTables->searchQuantsElectron[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);
                                                    
        // if(gridData->PressureE[i] < 1e-10) //Happens if the Table does not support the temperature
        // {
        //     gridData->PressureE[i] = 1e-10; //Hard limit 
        // }

        gridData->DpDtE[i] = (eq12 - eq11) / (std::get<1>(mTables->searchQuantsElectron[i].te_values) - std::get<0>(mTables->searchQuantsElectron[i].te_values));   
        gridData->DpDtE[i] = mTables->getGradientQuantity(std::get<0>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<0>(mTables->searchQuantsElectron[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].rho_values),
                                                    eq11, eq12, eq21, eq22, 
                                                    gridData->Density[i]);   

        if((gridData->DpDtE[i] < 0) &&(gridData->PressureE[i] < 0)) //Happens if the Table does not support the temperature
        {
            gridData->DpDtE[i] = 1e-10; //Hard limit 
        }
        gridData->PressureI[i] = mTables->biLinearInterpolation(std::get<0>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<0>(mTables->searchQuantsIon[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].rho_values),
                                                    iq11, iq12, iq21, iq22, gridData->TemperatureI[i], 
                                                    gridData->Density[i]);

        // if(gridData->PressureI[i] < 1e-10) //Happens if the Table does not support the temperature
        // {
        //     gridData->PressureI[i] = 1e-10; //Hard limit 
        // }
        // gridData->DpDtI[i] = (iq12 - iq11) / (std::get<1>(mTables->searchQuantsIon[i].te_values) - std::get<0>(mTables->searchQuantsIon[i].te_values));   
        gridData->DpDtI[i] = mTables->getGradientQuantity(std::get<0>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<0>(mTables->searchQuantsIon[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].rho_values),
                                                    iq11, iq12, iq21, iq22, 
                                                    gridData->Density[i]);

        if((gridData->DpDtI[i] < 0) &&(gridData->PressureI[i] < 0)) //Happens if the Table does not support the temperature
        {
            gridData->DpDtI[i] = 1e-10; //Hard limit 
        }
        gridData->TotalPressure[i] = gridData->PressureE[i] + gridData->PressureI[i];
        if(gridData->TotalPressure[i] < mLocalPressureCap)
        {
            gridData->TotalPressure[i] = mLocalPressureCap;
        }

    }

}
void LoadInEoS::updateIntEnergyTerms(GridData *gridData)
{
    for(int i = mStartNx; i < mNx; i++)
    {

        int ncols = mTables->FEOSTemperature.size();
        int elower_temp_index = std::get<0>(mTables->searchQuantsElectron[i].te_index);
        int eupper_temp_index = std::get<1>(mTables->searchQuantsElectron[i].te_index);

        int elower_density_index = std::get<0>(mTables->searchQuantsElectron[i].rho_index);
        int eupper_density_index = std::get<1>(mTables->searchQuantsElectron[i].rho_index);
        
        double eq11 = mIntMultiE*mTables->intEnergyElectronFeosTable[elower_density_index * ncols + elower_temp_index]; //q11
        double eq12 = mIntMultiE*mTables->intEnergyElectronFeosTable[elower_density_index * ncols + eupper_temp_index]; //q12
        double eq21 = mIntMultiE*mTables->intEnergyElectronFeosTable[eupper_density_index * ncols + elower_temp_index]; //q21
        double eq22 = mIntMultiE*mTables->intEnergyElectronFeosTable[eupper_density_index * ncols + eupper_temp_index]; //q22

        int ilower_temp_index = std::get<0>(mTables->searchQuantsIon[i].te_index);
        int iupper_temp_index = std::get<1>(mTables->searchQuantsIon[i].te_index);

        int ilower_density_index = std::get<0>(mTables->searchQuantsIon[i].rho_index);
        int iupper_density_index = std::get<1>(mTables->searchQuantsIon[i].rho_index);
        
        double iq11 = mIntMultiI*mTables->intEnergyIonFeosTable[ilower_density_index * ncols + ilower_temp_index]; //q11
        double iq12 = mIntMultiI*mTables->intEnergyIonFeosTable[ilower_density_index * ncols + iupper_temp_index]; //q12
        double iq21 = mIntMultiI*mTables->intEnergyIonFeosTable[iupper_density_index * ncols + ilower_temp_index]; //q21
        double iq22 = mIntMultiI*mTables->intEnergyIonFeosTable[iupper_density_index * ncols + iupper_temp_index]; //q22
        
        
        gridData->InternalEnergyE[i] = mTables->biLinearInterpolation(std::get<0>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<0>(mTables->searchQuantsElectron[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].rho_values),
                                                    eq11, eq12, eq21, eq22, gridData->TemperatureE[i], 
                                                    gridData->Density[i]);
        
        ////CHANGE THIS OBVIOUSLY NOT RIGHT, if speciifc heat goes negative possible if table is not fully resolved go to Ideal GAS
        if(((eq12 - eq11) / (std::get<1>(mTables->searchQuantsElectron[i].te_values) - std::get<0>(mTables->searchQuantsElectron[i].te_values))) <= 0)
        {
            gridData->SpecificHeatE[i] = (gridData->NumberDensityE[i] * 1.383E-23) / (gridData->Density[i] * (1.6667 - 1)); //FIXME IMPLEMENT CORRECTLY
        }
        else
        {
            // gridData->SpecificHeatE[i] = (eq12 - eq11) / (std::get<1>(mTables->searchQuantsElectron[i].te_values) - std::get<0>(mTables->searchQuantsElectron[i].te_values));   
            gridData->SpecificHeatE[i] =  mTables->getGradientQuantity(std::get<0>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].te_values),
                                                    std::get<0>(mTables->searchQuantsElectron[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsElectron[i].rho_values),
                                                    eq11, eq12, eq21, eq22, 
                                                    gridData->Density[i]);
        }
        
        

        gridData->InternalEnergyI[i] =  mTables->biLinearInterpolation(std::get<0>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<0>(mTables->searchQuantsIon[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].rho_values),
                                                    iq11, iq12, iq21, iq22, gridData->TemperatureI[i], 
                                                    gridData->Density[i]);

        if(((iq12 - iq11) / (std::get<1>(mTables->searchQuantsIon[i].te_values) - std::get<0>(mTables->searchQuantsIon[i].te_values))) <= 0)
        {
            gridData->SpecificHeatI[i] = (gridData->NumberDensityI[i] * 1.383E-23) / (gridData->Density[i] * (1.6667 - 1)); //FIXME IMPLEMENT CORRECTLY
        }
        else
        {
            gridData->SpecificHeatI[i] =  mTables->getGradientQuantity(std::get<0>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].te_values),
                                                    std::get<0>(mTables->searchQuantsIon[i].rho_values),
                                                    std::get<1>(mTables->searchQuantsIon[i].rho_values),
                                                    iq11, iq12, iq21, iq22, 
                                                    gridData->Density[i]);
        }
    }
}

void LoadInEoS::setMultipliers(float sphecamE, float sphecamI, float premE, float premI)
{
    
    mIntMultiE = sphecamE;
    mIntMultiI = sphecamI;
    mPreMultiE = premE;
    mPreMultiI = premI;
}


