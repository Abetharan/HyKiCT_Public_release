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
#ifndef INITHEADER
#define INITHEADER
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include "ExtSource.h"
#include "IntSource.h"
#include "FluidDynamics.h"
#include "GridData.h"
#include "Switches.h"
#include "FixedData.h"
#include "IO.h"
#include "IdealGasEoS.h"
#include "LoadInEoS.h"
#include "Ionization.h"
#include "EoS.hpp"
#include "PlasmaParameters.h"
#include "TimeStep.h"
#include "RadTrans.hpp"
#include "MultiRadTrans.h"
#include "LoadInTables.h"
// #include "MatrixSolver.h"
#include "SNBSource.h"
class Init
{             
    private:
        int mfindChange(std::vector<double> Ar);
    public:
        std::vector<std::unique_ptr<EoS>> initEoS(FixedData const &fixedData, Switches const &switchData, IO const &ioVector);
        std::vector<Ionization> initIonization(FixedData const&fixedData, Switches const&switchData, IO const &ioData);
        std::vector<LoadInTables> initOpacityTables(FixedData const&fixedData, Switches const&switchData, IO const &ioData);
        std::vector<std::unique_ptr<RadTrans>> initRadiation(FixedData const&fixedData);
        void initRadEnergy(GridData *gridData, FixedData const &fixedData,std::vector<MultiRadTrans>  &radiationTransport);
        
        void initialise(GridData *gridData, FixedData const &fixedData, Switches const &switchData, std::vector<std::unique_ptr<EoS>> &EquationsOfState,
                    std::vector<Ionization> &IonizationModels, PlasmaParameters &plasmaParams, ExtSource &ExternalSources, 
                    IntSource &InternalSources,std::vector<SNBSource> &snbSource, FluidDynamics &Hydro, TimeStep const &timeData,std::vector<MultiRadTrans>  &radiationTransport, bool operator_split = false);

};
#endif 