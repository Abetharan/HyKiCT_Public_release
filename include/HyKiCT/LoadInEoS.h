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
#ifndef LOADINFEOS_HPP
#define LOADINFEOS_HPP
#include "EoS.hpp"
#include "GridData.h"
#include "FixedData.h"
#include "LoadInTables.h"
#include <utility>
#include <functional> 
class LoadInEoS : public EoS
{
        private:
            int mNx;
            int mStartNx;
            std::unique_ptr<LoadInTables> mTables;
            double mLocalPressureCap = 1e-10;
            double mLocalIntEnergyCap = 1e-20;
            float mIntMultiE, mIntMultiI, mPreMultiE, mPreMultiI;

        public:
            LoadInEoS();
            LoadInEoS(int startNx, int nx, std::unique_ptr<LoadInTables> &tables);
            // LoadInEoS(int startNx, int nx, std::string path);
            void setSearchIndicies(std::unique_ptr<LoadInTables> &table) override;
            void updatePressureTerms(GridData *gridData) override;
            void updateIntEnergyTerms(GridData *gridData) override;
            void FindIndicies(GridData *gridData) override;
            // void Multipliers(GridData *gridData, float sphecamE, float sphecamI, float premE, float premI) override;
            void setMultipliers(float sphecamE, float sphecamI, float premE, float premI) override;
};
#endif