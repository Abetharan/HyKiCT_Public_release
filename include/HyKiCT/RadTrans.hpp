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
#ifndef RADTRANS_HPP
#define RADTRANS_HPP
#include "GridData.h"
#include "FixedData.h"
#include "TimeStep.h"
#include "Switches.h"
// #include "MatrixSolver.h"
#include "LoadInTables.h"
class RadTrans
{
    public:
        virtual ~RadTrans() = default;
        virtual void updateRadiationQuantities(GridData const *gridDataT, FixedData const &fixedData, TimeStep const &timeData, int updateType) = 0;
        virtual void setMaterialInterface(int index) = 0;
        virtual void storeOpacityTables(std::vector<LoadInTables> table) = 0;
        virtual void setEnergyGroup(int group){std::cerr<<"Not implemented" << "\n";};
        virtual void gather(GridData *gridData){std::cerr<<"Not implemented" << "\n";};
};
#endif