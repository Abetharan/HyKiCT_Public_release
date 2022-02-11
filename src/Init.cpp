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
#include "HyKiCT/Init.h"

int Init::mfindChange(std::vector<double> Ar)
{
    for(unsigned long int i = 1; i < Ar.size(); i++)
    {
        if(Ar[i] != Ar[i - 1])
        {
            return i;
        }
        else
        {
            continue;
        }
    }
    return 0;
}

std::vector<std::unique_ptr<EoS>> Init::initEoS(FixedData const &fixedData, Switches const &switchData, IO const &ioData)
{
    std::vector<std::unique_ptr <EoS>> EquationsOfState;
    if(switchData.MultiMaterial)
    {
        int material_1_start = 0;
        int material_1_end = fixedData.materialInterface;// mfindChange(fixedData.Ar);
        std::cout << "Material Interface is at " << material_1_end << std::endl;
        int material_2_start = fixedData.materialInterface; //material_1_end
        int material_2_end = fixedData.Nx;
        if(switchData.IdealGas)
        {
            EquationsOfState.emplace_back(std::make_unique<IdealGasEoS>(IdealGasEoS(material_1_start, material_1_end, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT)));
            EquationsOfState.emplace_back(std::make_unique<IdealGasEoS>(IdealGasEoS(material_2_start, material_2_end, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT)));
        }
        else
        {
            // LoadInTables material1();
            // LoadInTables material2();
            std::unique_ptr<LoadInTables> material_1 = std::make_unique<LoadInTables>(LoadInTables(material_1_start, material_1_end)); 
            std::unique_ptr<LoadInTables> material_2 = std::make_unique<LoadInTables>(LoadInTables(material_2_start, material_2_end)); 
            material_1->feosEoSLoadIn(ioData.FEOSPathMaterial1);
            material_1->converTemperatureToeV(fixedData.ev_to_k);
            material_2->feosEoSLoadIn(ioData.FEOSPathMaterial2);
            material_2->converTemperatureToeV(fixedData.ev_to_k);
            EquationsOfState.emplace_back(std::make_unique<LoadInEoS>(LoadInEoS(material_1_start, material_1_end, material_1)));
            EquationsOfState.at(0)->setMultipliers(fixedData.sphecamE1, fixedData.sphecamI1, fixedData.prrm1, fixedData.prrm1);
            EquationsOfState.emplace_back(std::make_unique<LoadInEoS>(LoadInEoS(material_2_start, material_2_end, material_2)));
            EquationsOfState.at(1)->setMultipliers(fixedData.sphecamE2, fixedData.sphecamI2, fixedData.prrm2, fixedData.prrm2);
        }
    }
    else
    {
        if(switchData.IdealGas)
        {
            EquationsOfState.emplace_back(std::make_unique<IdealGasEoS>(IdealGasEoS(0, fixedData.Nx, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT)));
        }
        else
        {
            std::unique_ptr<LoadInTables> material_1 = std::make_unique<LoadInTables>(LoadInTables(0, fixedData.Nx)); 
            material_1->feosEoSLoadIn(ioData.FEOSPathMaterial1);
            material_1->converTemperatureToeV(fixedData.ev_to_k);
            EquationsOfState.emplace_back(std::make_unique<LoadInEoS>(LoadInEoS(0, fixedData.Nx, material_1)));
            EquationsOfState.at(0)->setMultipliers(fixedData.sphecamE1, fixedData.sphecamI1, fixedData.prrm1, fixedData.prrm1);
        }
    }  
    return(EquationsOfState);
}

std::vector<Ionization> Init::initIonization(FixedData const&fixedData, Switches const&switchData, IO const &ioData)
{
    std::vector<Ionization> ionizationModels;
    if(switchData.MultiMaterial)
    {
        int material_1_start = 0;
        int material_1_end = fixedData.materialInterface;// mfindChange(fixedData.Ar);
        std::cout << "Material Interface is at " << material_1_end << std::endl;
        int material_2_start = fixedData.materialInterface; //material_1_end
        int material_2_end = fixedData.Nx;
        if(switchData.IonizationTables)
        {
            std::unique_ptr<LoadInTables> material_1 = std::make_unique<LoadInTables>(LoadInTables(material_1_start, material_1_end)); 
            std::unique_ptr<LoadInTables> material_2 = std::make_unique<LoadInTables>(LoadInTables(material_2_start, material_2_end)); 
            material_1->feosIonizationLoadIn(ioData.FEOSPathMaterial1);
            material_1->converTemperatureToeV(fixedData.ev_to_k);
            material_2->feosIonizationLoadIn(ioData.FEOSPathMaterial2);
            material_2->converTemperatureToeV(fixedData.ev_to_k);
            ionizationModels.push_back(Ionization(material_1_start, material_1_end, material_1));
            ionizationModels.push_back(Ionization(material_2_start, material_2_end, material_2));
        }
        else
        {
            ionizationModels.push_back(Ionization(material_1_start, material_1_end));
            ionizationModels.push_back(Ionization(material_2_start, material_2_end));
        }
    }
    else
    {
        if(switchData.IonizationTables)
        {
            std::unique_ptr<LoadInTables> material_1 = std::make_unique<LoadInTables>(LoadInTables(0, fixedData.Nx)); 
            material_1->feosIonizationLoadIn(ioData.FEOSPathMaterial1);
            material_1->converTemperatureToeV(fixedData.ev_to_k);
            ionizationModels.push_back(Ionization(0, fixedData.Nx, material_1));
        }
        else
        {
            ionizationModels.push_back(Ionization(0, fixedData.Nx));
        }
    }
    return(ionizationModels); 
}
std::vector<LoadInTables> Init::initOpacityTables(FixedData const&fixedData, Switches const&switchData, IO const &ioData)
{

    std::vector<LoadInTables> opacityTables;
    if(switchData.MultiMaterial)
    {
        int material_1_start = 0;
        int material_1_end = fixedData.materialInterface;// mfindChange(fixedData.Ar);
        std::cout << "Material Interface is at " << material_1_end << std::endl;
        int material_2_start = fixedData.materialInterface; //material_1_end
        int material_2_end = fixedData.Nx;
        LoadInTables material_1(material_1_start, material_1_end); 
        LoadInTables material_2(material_2_start, material_2_end); 
        material_1.OpacityLoadIn(ioData.OpacityPathMaterial1);
        // material_1.converTemperatureToeV(fixedData.ev_to_k);
        material_2.OpacityLoadIn(ioData.OpacityPathMaterial2);
        if(fixedData.RadNg > 1)
        {
            material_1.PhotonGridLoadIn(ioData.OpacityPathMaterial1);
            material_2.PhotonGridLoadIn(ioData.OpacityPathMaterial2);
        }
        else
        {
            material_1.PhotonGrid = {1e-50, 1e50};
            material_2.PhotonGrid = {1e-50, 1e50};
        }
        
        // material_2.converTemperatureToeV(fixedData.ev_to_k);
        opacityTables.push_back(material_1);
        opacityTables.push_back(material_2);
    }
    else
    {
        LoadInTables material_1(0, fixedData.Nx); 
        material_1.OpacityLoadIn(ioData.OpacityPathMaterial1);
        if(fixedData.RadNg > 1)
        {
            material_1.PhotonGridLoadIn(ioData.OpacityPathMaterial1);
        }
        else
        {
             material_1.PhotonGrid = {1e-50, 1e50};
        }
        material_1.converTemperatureToeV(fixedData.ev_to_k);
        opacityTables.push_back(material_1);
    }
    return(opacityTables); 
}
void Init::initRadEnergy(GridData *gridData, FixedData const &fixedData, std::vector<MultiRadTrans>  &radiationTransport)
{
    int j = 0;
    for(auto &i:radiationTransport)
    {
        std::vector<double> tmp(gridData->TmpALlRadEnergy.begin(), gridData->TmpALlRadEnergy.begin() + fixedData.Nx); 
        i.setRadEnergyDensity(tmp);
        gridData->TmpALlRadEnergy.erase(gridData->TmpALlRadEnergy.begin(), gridData->TmpALlRadEnergy.begin() + fixedData.Nx);
    }
}
void Init::initialise(GridData *gridData, FixedData const &fixedData, Switches const &switchData, std::vector<std::unique_ptr<EoS>> &EquationsOfState,
                    std::vector<Ionization> &IonizationModels,PlasmaParameters &plasmaParams, ExtSource &ExternalSources,
                     IntSource &InternalSources, std::vector<SNBSource> &snbSource, FluidDynamics &Hydro, TimeStep const &timeData, std::vector<MultiRadTrans>  &radiationTransport, bool operator_split)
{
    Hydro.updateNumberDensityI(gridData, fixedData);
    Hydro.updateCellCentreCoords(gridData);
    Hydro.updateDxs(gridData);
    if(switchData.MultiMaterial)
    {
        if(switchData.IonizationTables)
        {
            IonizationModels.at(0).updateIonization(gridData, fixedData, switchData);
            IonizationModels.at(1).updateIonization(gridData, fixedData, switchData);
        }
        else
        {
            IonizationModels.at(0).updateIonization(gridData, fixedData, switchData);
            IonizationModels.at(1).updateIonization(gridData, fixedData, switchData);
        }
        
    }
    else
    {
        if(switchData.IonizationTables)
        {
            IonizationModels.at(0).updateIonization(gridData, fixedData, switchData);
        }
        else
        {
            IonizationModels.at(0).updateIonization(gridData, fixedData, switchData);
        }
    }

    Hydro.updateNumberDensityE(gridData);
    if(switchData.MultiMaterial)
    {
        if(switchData.IdealGas)
        {
            EquationsOfState.at(0)->updatePressureTerms(gridData);
            EquationsOfState.at(0)->updateIntEnergyTerms(gridData);
            EquationsOfState.at(1)->updatePressureTerms(gridData);
            EquationsOfState.at(1)->updateIntEnergyTerms(gridData);
        }
        else
        {
            if(!switchData.IonizationTables)
            {
                EquationsOfState.at(0)->FindIndicies(gridData);
                EquationsOfState.at(1)->FindIndicies(gridData);
            }
            else
            {
                EquationsOfState.at(0)->setSearchIndicies(IonizationModels.at(0).mIonizationTable);
                EquationsOfState.at(1)->setSearchIndicies(IonizationModels.at(1).mIonizationTable);
            }
            
            EquationsOfState.at(0)->updatePressureTerms(gridData);
            EquationsOfState.at(0)->updateIntEnergyTerms(gridData);
            EquationsOfState.at(1)->updatePressureTerms(gridData);
            EquationsOfState.at(1)->updateIntEnergyTerms(gridData);
        }
        
    }    
    else
    {
        if(switchData.IdealGas)
        {
            EquationsOfState.at(0)->updatePressureTerms(gridData);
            EquationsOfState.at(0)->updateIntEnergyTerms(gridData);
        }
        else
        {
            if(!switchData.IonizationTables)
            {
                EquationsOfState.at(0)->FindIndicies(gridData);
            }
            else
            {
                EquationsOfState.at(0)->setSearchIndicies(IonizationModels.at(0).mIonizationTable);
            }
            
            EquationsOfState.at(0)->updatePressureTerms(gridData);
            EquationsOfState.at(0)->updateIntEnergyTerms(gridData);
        }
    }
    
    if(switchData.ConstantCoulombLog)
    {
        plasmaParams.setCoulombLog(gridData, 10);
    }
    if((!switchData.LeeMore) || (switchData.SNBHeatFlow))
    {
        plasmaParams.calculateCoulombLogEE(gridData, fixedData);
        plasmaParams.calculateCoulombLogEI(gridData->CoulombLogEI, gridData->TemperatureE,
                                    gridData->NumberDensityE, gridData->Zbar, fixedData);   
    }
    else
    {
        plasmaParams.calculateMinImpactParameter(gridData, fixedData);
        plasmaParams.calculateLeeMoreCoulombLogs(gridData, fixedData);
    }
    if(!operator_split)
    {
        if(switchData.InvBremsstrahlungOn)
        {
            if(switchData.LoadInLaserProfile)
            {
                ExternalSources.updateLaserPower(timeData.TotalTime);
            }
            // gridData->LinearInterpolate();
            // plasmaParams.calculateCoulombLogEI(gridData->InterCoulombLog, gridData->InterTemperatureE,
            //                             gridData->InterNumberDensityE, gridData->InterZbar, fixedData);   
            // ExternalSources.inverseBrem(gridData, fixedData);
            plasmaParams.calculatePlasmaFrequency(gridData, fixedData);
            plasmaParams.calculateMinImpactParameter(gridData, fixedData);
            plasmaParams.calculateThermalVelocity(gridData, fixedData);
            plasmaParams.calculateCoulombLogLaser(gridData, fixedData, ExternalSources.LaserWavelength);
            plasmaParams.calculateCollisionFrequencyEIOverC(gridData, fixedData, true);
            ExternalSources.calculateBeta(gridData);
            ExternalSources.inverseBrem(gridData);
        }
        if(switchData.HeatConductionOn)
        {
            InternalSources.calculateKappa(gridData, fixedData);
            if(!switchData.CoupleDivQ)
            {
                InternalSources.heatFlowE(gridData, fixedData, switchData);                    
                if(switchData.SNBHeatFlow)
                {
                    gridData->LinearInterpolate();
                    plasmaParams.calculateCoulombLogEI(gridData->InterCoulombLog, gridData->InterTemperatureE,
                                                gridData->InterNumberDensityE, gridData->InterZbar, fixedData);   
                    for(auto &i:snbSource)
                    {
                        i.snbHeatFlowCorrection(gridData, fixedData);
                        i.SNBGather(gridData);
                    }
                    snbSource.at(0).SNBCorrect(gridData);
                }
                if(switchData.CoupleOperatorSplit)
                {
                    gridData->CoupleOperatorSplitHeatFlowE = gridData->HeatFlowE; 
                }
                if(switchData.CoupleMulti)
                {
                    InternalSources.multiplierHeatFlowE(gridData, fixedData);                    
                }
                if(switchData.CoupleSubtract)
                {
                    InternalSources.subtractHeatFlowE(gridData);
                }
            }
            if(switchData.CoupleDivQ)
            {
                if(switchData.CoupleOperatorSplit)
                {
                    InternalSources.heatFlowE(gridData, fixedData, switchData);                    
                    gridData->CoupleOperatorSplitHeatFlowE = gridData->HeatFlowE; 
                    InternalSources.operatorSplitThermalConducE(gridData);
                }
                gridData->HeatConductionE = gridData->LoadInDivQ;;
            }
            else
            {
                if(switchData.CoupleOperatorSplit)
                {
                    InternalSources.operatorSplitThermalConducE(gridData);
                }
                InternalSources.thermalConducE(gridData);
            }
            InternalSources.heatFlowI(gridData, fixedData);
            InternalSources.thermalConducI(gridData);
        }
        if(switchData.ExchangeOn)
        {
            if((switchData.LeeMore) && (switchData.SNBHeatFlow))
            {
                plasmaParams.calculateMinImpactParameter(gridData, fixedData);
                plasmaParams.calculateLeeMoreCoulombLogs(gridData, fixedData);
            }
            plasmaParams.calculateCollisionFrequencyEI(gridData, fixedData, false);
            InternalSources.exchange(gridData, fixedData, timeData);        
        }
        if(switchData.RadTransportOn)
        {
            if(switchData.LoadInRadEnergy) 
            {
                initRadEnergy(gridData, fixedData, radiationTransport);
            }
            int j = 1;
            for(auto &i:radiationTransport)
            {
                if(switchData.LoadInRadEnergy)
                {
                    i.calculateOpacities(gridData, fixedData, false);
                    i.planckIntegral(gridData, fixedData);
                    i.calcFreeFreeEmission(gridData, fixedData);
                    i.calcFreeFreeAbsorption(gridData, fixedData);
                }
                else
                {
                    i.updateRadiationQuantities(gridData, fixedData, timeData);
                }
                i.gather(gridData);
                if((fixedData.Np == 1)&& (j < fixedData.RadNg))
                {
                    radiationTransport.at(j).copySearchEntries(i);
                }
                j++;
            }
            radiationTransport.at(0).updateRadiationTemperature(gridData, fixedData);
            radiationTransport.at(0).getOpacities(gridData);
        }
    }
}
