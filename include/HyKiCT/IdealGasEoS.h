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
#ifndef IDEALGASEOS_HPP
#define IDEALGASEOS_HPP
#include "EoS.hpp"
#include "GridData.h"
#include "FixedData.h"

class IdealGasEoS : public EoS
{
        private:
            int mNx;
            int mStartNx; 
            float mGamma;
            double mBoltzmannConstant;
        public:
            IdealGasEoS();
            IdealGasEoS(int startNx, int nx, double gamma, double kb);
            double SpecificHeatCapacity(double ne, double rho, double kb);
            void updatePressureTerms(GridData *gridData) override;
            void updateIntEnergyTerms(GridData *gridData) override;
            void updateSpecificHeatCapacity(GridData *gridData, int i);
            void updateSpecificHeatCapacity(GridData *gridData);
            void updateMarshackSpecificHeatCapacity(GridData *gridData, int i);
            // void updateDpDt(DataHandler &data);
};
#endif