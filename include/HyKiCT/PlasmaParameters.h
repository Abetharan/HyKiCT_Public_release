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
#ifndef PLASMAPARAMETERS_HPP
#define PLASMAPARAMETERS_HPP
#include <vector>
#include <algorithm>
#include "FixedData.h"
#include "GridData.h"
#include "FastMaths.hpp"

class PlasmaParameters
{   
    private:
        int mNx;
        std::vector<double> mCoulombLog{};
    public:
        PlasmaParameters(int nx);
        void calculateCoulombLogEI(GridData *gridData, FixedData const &fixedData);
        void calculateCoulombLogEI(std::vector<double> &CoulombLogEI, 
                                std::vector<double> const &TemperatureE, 
                                std::vector<double> const &NumberDensityE,
                                std::vector<double> const &Zbar, FixedData const &fixedData);
        void calculateCoulombLogEE(GridData *gridData, FixedData const &fixedData);
        void setCoulombLog(GridData *gridData, double coloumblog);
        void calculatePlasmaFrequency(GridData *gridData, FixedData const &fixedData);
        void calculateMinImpactParameter(GridData *gridData, FixedData const &fixedData);
        void calculateThermalVelocity(GridData *gridData, FixedData const &fixedData);
        void calculateCollisionFrequencyEI(GridData *gridData, FixedData const &fixedData, bool laser);
        void calculateCollisionFrequencyEIOverC(GridData *gridData, FixedData const &fixedData, bool laser);
        void calculateCoulombLogLaser(GridData *gridData, FixedData const &fixedData, double laserWavelength);
        void calculateLeeMoreCoulombLogs(GridData *gridData, FixedData const &fixedData);
};
#endif