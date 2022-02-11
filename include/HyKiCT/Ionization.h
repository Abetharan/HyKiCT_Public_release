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
#ifndef IONIZATION_HPP
#define IONIZATION_HPP
#include <memory>
#include "Switches.h"
#include "GridData.h"
#include "FixedData.h"
#include "LoadInTables.h"
class Ionization
{
    private:
        int mNx;
        int mStartNx;
        void mAproxThomasFermiIonize(GridData *gridData, FixedData const &fixedData);
        void mFullyIonized(GridData *gridData, FixedData const &fixedData);
        void mFEOSUpdateIonisation(GridData *gridData);
    public:
        Ionization(int startNx, int nx);
        Ionization(int startNx, int nx, std::unique_ptr<LoadInTables> &table);
        void updateIonization(GridData *gridData, FixedData const &fixedData,Switches const &switchData);
        std::unique_ptr<LoadInTables> mIonizationTable;
        std::unique_ptr<LoadInTables> getTable();

};
#endif