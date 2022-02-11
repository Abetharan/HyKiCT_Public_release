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
#include "gtest/gtest.h"
#include "HyKiCT/FixedData.h"
#include "HyKiCT/FluidDynamics.h"
#include "HyKiCT/GridData.h"
#include "HyKiCT/PlasmaParameters.h"
#include "HyKiCT/TimeStep.h"
#include "HyKiCT/Switches.h"
#include "HyKiCT/IdealGasEoS.h"
#include "HyKiCT/IO.h"


void initSodShockTube(GridData *gridData,int nx, float split, float gamma, double L, double mp, double kb)
{
    for(int i = 0; i < nx; i++)
    {
        if(i < nx*split )
        {
            gridData->TotalPressure[i] = 1;
            gridData->Density.push_back(1);
        }
        else
        {

            gridData->TotalPressure[i] = 0.1;
            gridData->Density.push_back(0.125);
        }
        gridData->Zbar[i] = 1;
        gridData->TemperatureE.push_back(0.0); // gridData->TotalPressure[i]/((gamma - 1 )* gridData->Density[i]));
        gridData->TemperatureI.push_back((mp * gridData->TotalPressure[i]) / (kb * gridData->Density[i]));
    }
    double dx = L / (nx + 1);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i*dx);
        gridData->Velocity.push_back(0);
    }
    for(int i = 0; i < nx; i++)
    {
        gridData->Mass.push_back(gridData->Density[i] * (gridData->CellWallCoord[i + 1] - gridData->CellWallCoord[i]));
    }
}
void loadSodExpected(std::string path, std::vector<double> &v,std::vector<double> &rho, std::vector<double> &pt, std::vector<double> &coord)
{
    std::vector<std::string> vars = {"velocity.txt", "rho.txt", "pt.txt", "coord.txt"};
    int no_var = vars.size();
    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::string file_path = path + "/" +  vars[i];
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    v.push_back(stod(strInput));
                }
                if(i == 1)
                {

                    rho.push_back(stod(strInput));
                }
                if(i == 2)
                {
                    pt.push_back(stod(strInput));
                }
                if(i == 3)
                {
                    coord.push_back(stod(strInput));
                }
            }
        }
    }
}
TEST(FluidDynamics, DISABLED_SodShockTube)
{
    int nx = 2000;
    float split = 0.5;
    float gamma = 1.66666667;
    double L = 1;
    TimeStep timeData(nx);
    timeData.Dt05 = 1e-7;
    timeData.Dt1 = 1e-7;
    timeData.Tmax = 0.2;
    timeData.dtGlobalMax = 1e-2;
    timeData.dtGlobalMin = 1e-7;
    timeData.setTimeParameters(0,0,0,0,0,0.1);
    FixedData fixedData;
    fixedData.Nx = nx;
    fixedData.RightFluidBoundaryCondition = "r";
    FluidDynamics Hydro(nx);
    IdealGasEoS EoS(0, nx, gamma, fixedData.BOLTZMANN_CONSTANT);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    std::vector<GridData> variableData;
    GridData initGridData(nx,true);
    initSodShockTube(&initGridData, nx, split, gamma, L, fixedData.PROTON_MASS, fixedData.BOLTZMANN_CONSTANT);
    fixedData.Z.resize(nx, 1);
    fixedData.Ar.resize(nx, 1);
    Hydro.updateDxs(&initGridData);
    Hydro.updateNumberDensityI(&initGridData,fixedData);
    Hydro.updateNumberDensityE(&initGridData);
    Hydro.updateCellCentreCoords(&initGridData);
    EoS.updateSpecificHeatCapacity(&initGridData);
    EoS.updateIntEnergyTerms(&initGridData);
    variableData.push_back(initGridData);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/SOD_SHOCKTUBE";
    std::vector<double> v, rho,pt, coord;
    loadSodExpected(localPath, v, rho, pt, coord);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {
            // variableData.back().dump(ioData, timeData.Step);
            break;
        }
        if(variableData.size() > 2)
        {
                variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        //Units being tested

        variableData.at(lastStepDataIndex).Mass = variableData.at(lastStepDataIndex - 1).Mass;
        variableData.at(lastStepDataIndex).Zbar = variableData.at(lastStepDataIndex - 1).Zbar;
        
        Hydro.updateViscosity(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), fixedData);
        Hydro.updateMomentum(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), fixedData, timeData);
        Hydro.updateCoords(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData);
        Hydro.updateCellCentreCoords(&variableData.at(lastStepDataIndex));
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateDensity(&variableData.at(lastStepDataIndex));

        Hydro.updateNumberDensityI(&variableData.at(lastStepDataIndex), fixedData);
        Hydro.updateNumberDensityE(&variableData.at(lastStepDataIndex));
        Hydro.updateTemperatureI(&variableData.at(lastStepDataIndex- 1), &variableData.at(lastStepDataIndex), timeData, true);
        variableData.at(lastStepDataIndex).TemperatureE.resize(nx, 1);
        EoS.updatePressureTerms(&variableData.at(lastStepDataIndex));
        EoS.updateSpecificHeatCapacity(&variableData.at(lastStepDataIndex));
        EoS.updateIntEnergyTerms(&variableData.at(lastStepDataIndex));

        // if(timeData.Step % 100 == 0)
        // {
        //     variableData.back().dump(ioData, timeData.Step);
        // }
        timeData.incrementStep();
        timeData.setTotalTime();
        timeData.calculateNewTime(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex));    
    }

    std::vector<double> hy_rho = variableData.back().Density; 
    std::vector<double> hy_vel = variableData.back().Velocity;
    std::vector<double> hy_pt = variableData.back().TotalPressure;
    std::vector<double> hy_x = variableData.back().CellWallCoord;
    std::vector<double> relative_errors_u0;
    std::vector<double> relative_errors_u1;
    std::vector<double> relative_errors_v;
    std::vector<double> relative_errors_x;
    for(int i = 0; i < fixedData.Nx; i++)
    {
        relative_errors_u0.push_back(abs(hy_rho[i] - rho[i])/rho[i]);
        relative_errors_u1.push_back(abs(hy_pt[i] - pt[i])/pt[i]);
        if(v[i] == 0)
        {
            relative_errors_v.push_back(0.0);
        }
        else
        {
            relative_errors_v.push_back(abs(hy_vel[i] - v[i])/v[i]);
        }
        if(coord[i] == 0)
        {
            relative_errors_x.push_back(0.0);
        }
        else
        {
            relative_errors_x.push_back(abs(hy_x[i] - coord[i])/coord[i]);
        }


    } 
    auto max_error1 = std::max_element(std::begin(relative_errors_u1), std::end(relative_errors_u1));
    auto max_error2 = std::max_element(std::begin(relative_errors_u0), std::end(relative_errors_u0));
    auto max_error3 = std::max_element(std::begin(relative_errors_v), std::end(relative_errors_v));
    auto max_error4 = std::max_element(std::begin(relative_errors_x), std::end(relative_errors_x));
    // ASSERT_LT(*max_error1, 1.0e-8);
    // ASSERT_LT(*max_error2, 1.0e-8);
    // ASSERT_LT(*max_error3, 1.0e-8);
    // ASSERT_LT(*max_error4, 1.0e-8);
}
