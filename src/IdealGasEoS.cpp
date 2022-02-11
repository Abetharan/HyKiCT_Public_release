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
#include "HyKiCT/IdealGasEoS.h"

IdealGasEoS::IdealGasEoS(int startNx, int nx, double gamma, double kb)
{    
   mStartNx = startNx;
   mNx = nx;
   mGamma = gamma;
   mBoltzmannConstant = kb;

}
void IdealGasEoS::updatePressureTerms(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->PressureE[i] = gridData->NumberDensityE[i] * mBoltzmannConstant * gridData->TemperatureE[i];
        gridData->PressureI[i] = gridData->NumberDensityI[i] * mBoltzmannConstant * gridData->TemperatureI[i];
        gridData->TotalPressure[i] = gridData->PressureE[i] + gridData->PressureI[i];
        gridData->DpDtE[i] = gridData->NumberDensityE[i] * mBoltzmannConstant;
        gridData->DpDtI[i] = gridData->NumberDensityI[i] * mBoltzmannConstant;
        // For debug
        // std::cout << gridData->TemperatureE [i]<< " ";
        // std::cout << gridData->NumberDensityE[i] << " ";
        // std::cout << gridData->PressureE[i] << " ";
        // std::cout << "i" <<"\n";
    }
}
double IdealGasEoS::SpecificHeatCapacity(double n, double rho, double kb)
{
    return ((n *kb) / (rho * (mGamma - 1)));
}

void IdealGasEoS::updateSpecificHeatCapacity(GridData *gridData, int i)
{
        gridData->SpecificHeatE[i] = IdealGasEoS::SpecificHeatCapacity(gridData->NumberDensityE[i], gridData->Density[i], mBoltzmannConstant);
        gridData->SpecificHeatI[i] = IdealGasEoS::SpecificHeatCapacity(gridData->NumberDensityI[i], gridData->Density[i], mBoltzmannConstant);
}
void IdealGasEoS::updateSpecificHeatCapacity(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->SpecificHeatE[i] = (gridData->NumberDensityE[i] * mBoltzmannConstant) / (gridData->Density[i] * (mGamma - 1));
        gridData->SpecificHeatI[i] = (gridData->NumberDensityI[i] * mBoltzmannConstant) / (gridData->Density[i] * (mGamma - 1));
    }
}
void IdealGasEoS::updateIntEnergyTerms(GridData *gridData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->InternalEnergyE[i] = (gridData->NumberDensityE[i] * mBoltzmannConstant * gridData->TemperatureE[i]) / 
                                                (gridData->Density[i] * (mGamma - 1));

        gridData->InternalEnergyI[i] = (gridData->NumberDensityI[i] * mBoltzmannConstant * gridData->TemperatureI[i]) / 
                                                (gridData->Density[i] * (mGamma - 1));

        updateSpecificHeatCapacity(gridData , i);
    }
}
