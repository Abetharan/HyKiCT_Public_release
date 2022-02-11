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
#ifndef UPDATEHYDROHEADERDEF
#define UPDATEHYDROHEADERDEF
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <assert.h>
#include "FixedData.h"
#include "GridData.h"
#include "Switches.h"
#include "TimeStep.h"
// #include "MatrixSolver.h"
class FluidDynamics
{
private:
    //Const Value
    int mNx;
    void mSolveTridiagonal(double* a, double* b, double* c, double* d, int n);
    double mVelocityCap = 1e-15, mViscosityCap = 1e-15;
public:
    FluidDynamics(int nx);
    void updateViscosity(GridData *gridDataT_1,  GridData *gridDataT, FixedData const &fixedData);
    void updateMomentum(GridData *gridDataT_1,  GridData *gridDataT, FixedData const &fixedData, TimeStep const &timeData);
    void updateDensity(GridData *gridDataT);
    void updateCoords(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData);
    void updateCellCentreCoords(GridData *gridDataT);
    void updateDxs(GridData *gridDataT);
    void updateNumberDensityE(GridData *gridDataT);
    void updateNumberDensityI(GridData *gridDataT, FixedData const &fixedData);
    void updateTemperatureE(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic);
    void updateTemperatureEOperatorSplit(GridData *gridDataT, TimeStep const &timeData, std::vector<double> *source);
    void updateTemperatureEOperatorSplit(GridData *gridDataT_1, GridData *gridDataT, TimeStep const &timeData);
    void updateTemperatureCoupleOperatorSplit(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic);

    void updateTemperatureImplicitHeatConduction(std::vector<double> &Temperature, std::vector<double> &HeatKappa, std::vector<double> const &SpecificHeat, std::vector<double> const &Mass,
                                                std::vector<double> const &Dx_k, double dt05, std::vector<double> const *multipliers = nullptr, bool CoupleMulti = false);
    // void updateTemperatureImplicitHeatConduction(GridData *gridDataT, TimeStep const &timeData, Switches const &switchData);
    // void updateTemperatureIImplicitHeatConduction(GridData *gridDataT, TimeStep const &timeData, Switches const &switchData);
    void updateTemperatureI(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic);
    void updateTemperatureIOperatorSplit(GridData *gridDataT, TimeStep const &timeData, std::vector<double> *source);
    void updateTemperatureIOperatorSplit(GridData *gridDataT_1, GridData *gridDataT, TimeStep const &timeData);
    void fluidLeakSources(GridData *gridData, FixedData const &fixedData, double area, double start_index, double end_index);
};
#endif