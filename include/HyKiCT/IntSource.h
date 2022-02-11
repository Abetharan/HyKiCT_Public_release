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
#ifndef INTSOURCE_HPP
#define INTSOURCE_HPP
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "GridData.h"
#include "FixedData.h"
#include "Switches.h"
#include "TimeStep.h"
#include "FastMaths.hpp"
class IntSource
{
    private:
        int mNx;
        double mCalculateKappaE(double TemperatureE, double CoulombLogEI, double Zbar);
        double mCalculateKappaI(double TemperatureI, double CoulombLogEI, double Zbar, int Ar);

    public:
        IntSource(int nx);
        void calculateKappa(GridData *gridData, FixedData const &fixedData);
        void calculateImplicitKappaE(GridData *gridData, FixedData const &fixedData);
        void calculateImplicitKappaI(GridData *gridData, FixedData const &fixedData);
        void calculateLeeMoreKappaECorr(GridData *gridData, FixedData const &fixedData);
        void heatFlowI(GridData *gridData, FixedData const &fixedData);
        void heatFlowE(GridData *gridData, FixedData const &fixedData, Switches const &switchData);
        void diagonsticHeatFlowE(GridData *gridData, FixedData const &fixedData, Switches const &switchData);
        void subtractHeatFlowE(GridData *gridData);
        void multiplierHeatFlowE(GridData *gridData, FixedData const &fixedData);
        void spitzerHarmHeatFlowE(GridData *gridData, FixedData const &fixedData);
        void thermalConducI(GridData *gridData);
        void thermalConducE(GridData *gridData);
        void operatorSplitThermalConducE(GridData *gridData);
        void exchange(GridData *gridData, FixedData const &fixedData, TimeStep const&timeData);
        void setExchangeCoefficient(GridData *gridData, double ex);

};



#endif