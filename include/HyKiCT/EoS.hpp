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
#ifndef EOSHEADERDEF
#define EOSHEADERDEF
#include <iostream>
#include "FixedData.h"
#include "GridData.h"
#include "LoadInTables.h"
class EoS
{
   public:
         virtual ~EoS() = default;
         virtual void updatePressureTerms(GridData *gridData)  = 0; //Update Pressure and dPdT
         virtual void updateIntEnergyTerms(GridData *gridData) = 0; //Update Int Energy and Cv(de/dT)
         virtual void setSearchIndicies(std::unique_ptr<LoadInTables> &table){std::cerr<<"Not implemented" << "\n";};
         virtual void FindIndicies(GridData *gridData) {std::cerr<<"Not implemented" << "\n";};
         virtual void setMultipliers(float sphecamE, float sphecamI, float premE, float premI){std::cerr<<"Not implemented" << "\n";};
          // virtual void updateDpDt() = 0;
          // virtual void updateSpecificHeatCapacity() = 0;
};
#endif