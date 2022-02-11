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
#include "HyKiCT/Switches.h"

Switches::Switches(std::string path)
{
    YAML::Node config = YAML::LoadFile(path);
    ViscosityOn = config["Switches"]["Viscosity"].as<bool>();
    VelocityOn = config["Switches"]["Velocity"].as<bool>();
    HeatConductionOn = config["Switches"]["HeatConduction"].as<bool>();
    ImplicitHeatConduction = config["Switches"]["ImplicitHeatConduction"].as<bool>();
    ExchangeOn = config["Switches"]["Exchange"].as<bool>();
    InvBremsstrahlungOn = config["Switches"]["InvBrem"].as<bool>();
    IsothermalModeON = config["Switches"]["IsothermalMode"].as<bool>();
    SingleTemperature = config["Switches"]["SingleTemperature"].as<bool>();
    MultiMaterial = config["Switches"]["MultiMaterial"].as<bool>();
    IdealGas = config["Switches"]["IdealGas"].as<bool>();
    EoSTables = config["Switches"]["EoSTables"].as<bool>();
    IonizationTables = config["Switches"]["IonizationTables"].as<bool>();
    AdapativeTimeStep = config["Switches"]["AdapativeTimeStep"].as<bool>();
    CoupleDivQ = config["Switches"]["CoupleDivQ"].as<bool>();
    CoupleMulti = config["Switches"]["CoupleMulti"].as<bool>();
    CoupleOperatorSplit = config["Switches"]["CoupleOperatorSplit"].as<bool>();
    CoupleSubtract= config["Switches"]["CoupleSubtract"].as<bool>();
    ConstantCoulombLog = config["Switches"]["ConstantCoulombLog"].as<bool>();
    RadTransportOn = config["Switches"]["RadiationTransport"].as<bool>();
    RadFFOn = config["Switches"]["RadFreeFree"].as<bool>();
    RadFicksLaw = config["Switches"]["RadDiffusion"].as<bool>();
    FullyIonized = config["Switches"]["FullyIonized"].as<bool>();
    LoadInLaserProfile = config["Switches"]["LoadInLaserProfile"].as<bool>();
    SNBHeatFlow = config["Switches"]["SNBHeatFlow"].as<bool>();
    LeeMore = config["Switches"]["LeeMore"].as<bool>();
    LoadInRadEnergy = config["Switches"]["LoadInRadEnergy"].as<bool>();
    LeakSource = config["Switches"]["LeakSource"].as<bool>();
    normalOrCouple();
    std::cout << "SWITCHES BEING USED" << "\n";
    std::cout << "\n";
    std::cout << "FLUID DYNAMIC SWITCHES:" << "\n";
    std::cout << "VISCOSITY: " << ViscosityOn << "\n";
    std::cout << "VELOCITY: " << VelocityOn << "\n";
    std::cout << "\n";
    std::cout << "EQUATION OF STATE SWITCHES:" << "\n";
    std::cout << "IDEAL GAS: " << IdealGas <<"\n";
    std::cout << "EOS Tables: " << EoSTables <<"\n";
    std::cout << "Ionization Tables: " << IonizationTables <<"\n";
    std::cout << "MULTI-MATERIAL: " << MultiMaterial <<"\n";
    std::cout << "\n";
    std::cout << "RADIATION TRANSPORT SWITCHES:" << "\n";
    std::cout << "RADIATION TRANSPORT: " <<RadTransportOn << "\n";
    std::cout << "RADIATION FREE-FREE EMISSION: " << RadFFOn<< "\n";
    std::cout << "RADIATION DIFFUSION: " << RadFicksLaw<<"\n";
    std::cout << "\n";
    std::cout << "COUPLED SWITCHES: " << "\n";
    std::cout << "COUPLING USING DIV Q: " <<CoupleDivQ << "\n";
    std::cout << "COUPLING USING MULTIPLIERS: "<< CoupleMulti << "\n";
    std::cout << "COUPLING VIA OPERATOR SPLIT: "<< CoupleOperatorSplit << "\n";
    std::cout << "\n";
    std::cout << "OTHER SWITCHES " << "\n";
    std::cout << "CONSTANT COULOMB LOG: " <<ConstantCoulombLog <<"\n";
    std::cout << "ADAPATIVE TIME STEPPING: " << AdapativeTimeStep<<"\n";
    std::cout << "SINGLE TEMPERATURE MODE: " <<SingleTemperature<<"\n";
    std::cout << "ISOTHERMAL MODE: " << IsothermalModeON<<"\n";
    std::cout << "ADIABATIC MODE: " << AdibaticModeON<<"\n";
    std::cout << std::endl;

} 
void Switches::normalOrCouple()
{
    if((CoupleDivQ) || (CoupleMulti) || (CoupleSubtract))
    {
        Couple = true;
    }
    else
    {
        Couple = false;
    }
    
}
