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
#ifndef SWITCHES_HPP
#define SWITCHES_HPP
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <yaml-cpp/yaml.h>

class Switches
{
    public:
        bool AdapativeTimeStep;
        bool AdibaticModeON = false; //REF obsolete remove
        bool ConstantCoulombLog;
        bool Couple; //not loaded in ... infered.
        bool CoupleDivQ;
        bool CoupleMulti;
        bool CoupleOperatorSplit;
        bool CoupleSubtract;
        bool ExchangeOn;
        bool FullyIonized;
        bool HeatConductionOn;
        bool SNBHeatFlow;
        bool LeeMore;
        bool ImplicitHeatConduction;
        bool LoadInRadEnergy;
        bool IdealGas;
        bool EoSTables;
        bool IonizationTables;
        bool InvBremsstrahlungOn;
        bool IsothermalModeON;
        bool LoadInLaserProfile;
        bool MultiMaterial;
        bool RadFFOn;
        bool RadFicksLaw;
        bool RadTransportOn;
        bool SingleTemperature;
        bool VelocityOn;
        bool ViscosityOn;
        bool LeakSource;
        Switches() {}; 
        Switches(std::string path);
        void normalOrCouple();
};

#endif
