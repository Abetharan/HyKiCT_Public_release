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
#define VERSION_MAJOR 2
#define VERSION_MINOR 4
#define VERSION_REVISION 7
#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <fenv.h>
#include <assert.h>
#include <algorithm>
// #include <petscksp.h>
#include "HyKiCT/IO.h"
#include "HyKiCT/GridData.h"
#include "HyKiCT/FixedData.h"
#include "HyKiCT/Switches.h"
#include "HyKiCT/Init.h"
#include "HyKiCT/FluidDynamics.h"
#include "HyKiCT/IntSource.h"
#include "HyKiCT/ExtSource.h"
#include "HyKiCT/EoS.hpp"
#include "HyKiCT/Ionization.h"
#include "HyKiCT/RadTrans.hpp"
#include "HyKiCT/MultiRadTrans.h"
#include "HyKiCT/TimeStep.h"
#include "HyKiCT/PlasmaParameters.h"
// #include "HyKiCT/MatrixSolver.h"
#include "HyKiCT/RadTrans.hpp"
// #include "HyKiCT/PhysicsThreads.hpp"
#include "HyKiCT/SNBSource.h"
static char help[] = "HyKiCT - Hydro-Kinetic Coupling Testbed"; 

static void showUsage(std::string name)
{
    std::cerr   << "Usage: " << name << " <option(s)> SOURCES"
                << " Options:\n"
                <<"\t-v,--version\tVersion\n"
                << "\t-h,--help\tShow this help message\n"
                << "\t-p,--path \tSpecify the Base Directory where Source Code folder is located\n  "
                << "\t-Vb,--verbose \tVerbose Ouput"
                << std::endl;
    exit(1);
}
int findChange(std::vector<double> Ar)
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

int main(int argc, char* argv[])
{
    std::string path;
    bool verbose_output = false;
    if(argc == 1)
    {
        std::cerr << "REQUIRES PATH" << std::endl;
        exit(1);
    }
    for(int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if((arg == "-v") || (arg == "--version"))
        {
            std::cout << "ELH-1 VERSION: " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << "\n";
            return(1);
        }
        if((arg == "-Vb") || (arg == "--verbose"))
        {
            verbose_output = true;
            ++i;
            arg = argv[i];
        }
        if((arg == "-h") || (arg == "--help"))
        {
            showUsage(argv[0]);
        }
        else if((arg == "-p") || (arg == "--path"))
        {
            path = argv[i + 1];
            std::cout << "args " << argv[i] << std::endl;
            std::cout << "args " << argv[i+1] << std::endl;
            break;
        }
        else
        {
            std::cerr << "INVALID ARGUEMENT" << std::endl;
            exit(1);
        }
        
    }

    //init All objects
    IO ioData(path);
    FixedData constantData;
    Switches switchData(path); 
    // if((switchData.RadTransportOn)||(switchData.SNBHeatFlow)||(switchData.ImplicitHeatConduction))
    // {
        //Initialise PETSC 
        // PetscMPIInt size;
        // PetscInitialize(&argc, &argv, (char*)0, help);
        // MPI_Comm_size(PETSC_COMM_WORLD, &size);
        // if(size != 1)SETERRQ(PETSC_COMM_WORLD, 1, "Only valid for one processor at the moment");
    // }
    constantData.init(path);
    constantData.initVectors(ioData.InitPath, switchData.CoupleMulti);
    TimeStep timeData(constantData.Nx, constantData.Gamma);
    timeData.init(path);

    ioData.createDirectories();
    //Create Data vector
    std::vector<GridData> variableData;
    GridData initGridData(constantData.Nx, true);
    int coupleMethod;
    if(switchData.CoupleMulti)
    {
       coupleMethod = 1; 
    }
    else if(switchData.CoupleDivQ)
    {
        coupleMethod = 2;
    }
    else
    {
        coupleMethod = 3;
    }
    initGridData.init(ioData.InitPath, switchData.LoadInRadEnergy, switchData.Couple, coupleMethod);
    
    Init init;
    FluidDynamics Hydro(constantData.Nx);
    // RadTrans radiationTransport(constantData.Nx);
    PlasmaParameters plasmaParams(constantData.Nx);
    IntSource internalSources(constantData.Nx);
    ExtSource externalSources(path);
    // PhysicsThreads allThreads;
    // std::vector<MatrixSolver> matrixSolver(1, MatrixSolver(constantData.Nx));
    std::vector<SNBSource> snbSources(constantData.Ng, constantData.Nx);
    constantData.Np = 1; //REF Threading IS DISABLED ... DOES NOT GIVE ANY SPEED UP WITH CURRENT SETUP where groups are < 20 
    // if(constantData.Np > 1)
    // {
    //     allThreads.resize(constantData.Np);
    // }
    // if((switchData.RadTransportOn) || (switchData.SNBHeatFlow))
    // {
        // for(auto &i:matrixSolver)
        // {
        //     i.createPETScObjects();
        // }
    // }
    if(switchData.LoadInLaserProfile)
    {
        externalSources.Init(ioData.LaserProfilePath);
    }
    if(switchData.MultiMaterial)
    {
        constantData.materialInterface = findChange(constantData.Ar);
    }
    if(switchData.LeakSource)
    {
        constantData.loadLeakAreas(ioData.LeakProfilePath);
    }
    std::vector<std::unique_ptr<EoS>> equationsOfState = init.initEoS(constantData, switchData, ioData);
    std::vector<Ionization> ionizationModels = init.initIonization(constantData, switchData, ioData);

    std::vector<MultiRadTrans> radiationTransport(constantData.RadNg, MultiRadTrans(constantData.Nx));
    std::vector<LoadInTables> opacityTable;
    if(switchData.RadTransportOn)
    {

        opacityTable = init.initOpacityTables(constantData, switchData, ioData);
        int j = 0;
        for(auto &i:radiationTransport)
        {
            i.storeOpacityTables(opacityTable);
            i.setEnergyGroup(j);
            i.setTotalEnergyGroup(constantData.RadNg);
            if(switchData.MultiMaterial)
            {
                i.setMaterialInterface(constantData.materialInterface);
            }
            if(switchData.LeakSource)
            {
                i.setLeakArea(constantData, 0.0);
            }
            j++;
            if(j > constantData.RadNg)
            {
                std::cerr << "Too many groups" << std::endl;
            }
        }
    }
    if(switchData.SNBHeatFlow)
    {
     
        int j = 0;
        double max_energy_group = *max_element(initGridData.TemperatureE.begin(),
                                    initGridData.TemperatureE.end())*constantData.MaxE;
        for(auto &i:snbSources)
        {
            i.InitSNB(constantData.Ng, j,max_energy_group);
            j++;
        }
    }

    //Calcualte all t=0 quantites self-consistently from density,temperature and material properties
    init.initialise(&initGridData, constantData, switchData, equationsOfState, ionizationModels, plasmaParams, 
                    externalSources, internalSources, snbSources, Hydro, timeData, radiationTransport);
    initGridData.dump(ioData, timeData.Step);
    initGridData.terminalPrint();
    variableData.push_back(initGridData);
    timeData.dump(ioData); 

    //Set Initial DT
    int i = 1;
    while(true)
    {   
        //if t=0 init kill program immediatly after init
        if((timeData.MaxSteps== 0) && timeData.Tmax == 0)
        {
            return EXIT_SUCCESS;
        }
        //House Keeping Printing
        if((timeData.MaxSteps > 1))
        {
            //Last step dump
            if(timeData.MaxSteps - i  == 0)
            {
                std::cout << "Current Step is :" << i << " timeData Remaining : " << timeData.Tmax - timeData.TotalTime << ". With dt05 : " << timeData.Dt05 << " & dt01 : " << timeData.Dt1 << " Total timeData in seconds : " 
                        << timeData.TotalTime <<". \n";
                variableData.back().dump(ioData,timeData.Step);
                timeData.dump(ioData);
                break;

            }
        
            if((timeData.Step % constantData.OutputFrequency == 0))
            {
            std::cout << "Current Step is :" << i << " Steps Remaining : " << timeData.MaxSteps - i << ". With dt05 : " << timeData.Dt05 << " & dt01 : " << timeData.Dt1 << " Total timeData in seconds : " 
                        << timeData.TotalTime <<". \n";
            }
        }
        else
        {
            //last step dump
            if(timeData.TotalTime >= timeData.Tmax)
            {
                std::cout << "Current Step is :" << i << " timeData Remaining : " << timeData.Tmax - timeData.TotalTime << ". With dt05 : " << timeData.Dt05 << " & dt01 : " << timeData.Dt1 << " Total timeData in seconds : " 
                        << timeData.TotalTime <<". \n";
                variableData.back().dump(ioData, timeData.Step);
                timeData.dump(ioData);
                break;
            }
            if((timeData.Step % constantData.OutputFrequency == 0))
            {
                std::cout << "Current Step is :" << i << " timeData Remaining : " << timeData.Tmax - timeData.TotalTime << ". With dt05 : " << timeData.Dt05 << " & dt01 : " << timeData.Dt1 << " Total timeData in seconds : " 
                        << timeData.TotalTime <<". \n";
            }
        }
        
        //To reduce memory overhead.. keep data structor to be only last 4 time steps 
        if(variableData.size() > 3)
        {
                variableData.erase(variableData.begin());   
        } 
        //Init next gridData
        GridData grid(constantData.Nx);
        variableData.push_back(grid);
        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size

        variableData.at(lastStepDataIndex).Mass =  variableData.at(lastStepDataIndex - 1).Mass; //Mass stays fixed between steps 

        //Update Coulomb Logs
        
        //Hydro Step  
        if(switchData.ViscosityOn)
        {
            Hydro.updateViscosity(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), constantData);
        }
        if(switchData.VelocityOn)
        {
            Hydro.updateMomentum(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), constantData, timeData);
            Hydro.updateCoords(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData);
        }
        else
        {
            variableData.at(lastStepDataIndex).Velocity = variableData.at(lastStepDataIndex - 1).Velocity;
            variableData.at(lastStepDataIndex).CellWallCoord = variableData.at(lastStepDataIndex - 1).CellWallCoord;
        }
        
        Hydro.updateCellCentreCoords(&variableData.at(lastStepDataIndex));
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateDensity(&variableData.at(lastStepDataIndex));
        Hydro.updateNumberDensityI(&variableData.at(lastStepDataIndex), constantData);
        //If IsoThermal Mode is ON. DO NOT UPDATE TEMPERATURE! i.e. if is
        if(!switchData.IsothermalModeON)
        {
            //Electron
            if(switchData.CoupleOperatorSplit)
            {
                Hydro.updateTemperatureCoupleOperatorSplit(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, switchData.AdibaticModeON);
            }
            else
            {
                Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, switchData.AdibaticModeON);
            }

            if(switchData.SingleTemperature)
            {
                //electron temperature == ion temperature
                variableData.at(lastStepDataIndex).TemperatureI = variableData.at(lastStepDataIndex).TemperatureE;
            }
            else
            {   
                Hydro.updateTemperatureI(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, switchData.AdibaticModeON);
            }
        }
        else
        { //Constant temperature
            variableData.at(lastStepDataIndex).TemperatureE = variableData.at(lastStepDataIndex - 1).TemperatureE;
            variableData.at(lastStepDataIndex).TemperatureI = variableData.at(lastStepDataIndex - 1).TemperatureI;
        }

        //EoS Step
        if(switchData.MultiMaterial)
        {
            if(switchData.IonizationTables)
            {
                ionizationModels.at(0).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
                ionizationModels.at(1).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
            }
            else
            {
                ionizationModels.at(0).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
                ionizationModels.at(1).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
            }
            
        }
        else
        {
            if(switchData.IonizationTables)
            {
                ionizationModels.at(0).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
            }
            else
            {
                ionizationModels.at(0).updateIonization(&variableData.at(lastStepDataIndex), constantData, switchData);
            }
        }

        //Consisten update in ne AFTER calculating ionization.
        Hydro.updateNumberDensityE(&variableData.at(lastStepDataIndex));

    if(switchData.MultiMaterial)
    {
        if(switchData.IdealGas)
        {
            equationsOfState.at(0)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(0)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(1)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(1)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
        }
        else
        {
            if(!switchData.IonizationTables)
            {
                equationsOfState.at(0)->FindIndicies(&variableData.at(lastStepDataIndex));
                equationsOfState.at(1)->FindIndicies(&variableData.at(lastStepDataIndex));
            }
            else
            {
                equationsOfState.at(0)->setSearchIndicies(ionizationModels.at(0).mIonizationTable);
                equationsOfState.at(1)->setSearchIndicies(ionizationModels.at(1).mIonizationTable);
            }
            equationsOfState.at(0)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(0)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(1)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(1)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
        }
        
    }    
    else
    {
        if(switchData.IdealGas)
        {
            equationsOfState.at(0)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(0)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
        }
        else
        {
            if(!switchData.IonizationTables)
            {
                equationsOfState.at(0)->FindIndicies(&variableData.at(lastStepDataIndex));
            }
            else
            {
                equationsOfState.at(0)->setSearchIndicies(ionizationModels.at(0).mIonizationTable);
            }
            equationsOfState.at(0)->updatePressureTerms(&variableData.at(lastStepDataIndex));
            equationsOfState.at(0)->updateIntEnergyTerms(&variableData.at(lastStepDataIndex));
        }
    }
        //Calculate Plasma properties 
        //Calcualte coulomb logs
        if(switchData.ConstantCoulombLog)
        {
            plasmaParams.setCoulombLog(&variableData.at(lastStepDataIndex), 10);
        }
        else
        {
            plasmaParams.calculateCoulombLogEE(&variableData.at(lastStepDataIndex), constantData);
            plasmaParams.calculateCoulombLogEI(variableData.at(lastStepDataIndex).CoulombLogEI,
                                            variableData.at(lastStepDataIndex).TemperatureE,
                                            variableData.at(lastStepDataIndex).NumberDensityE, 
                                            variableData.at(lastStepDataIndex).Zbar,constantData);
        }

        //Laser heating
        if(switchData.InvBremsstrahlungOn)
        {
            if(timeData.TotalTime <= externalSources.LaserDuration)
            {
                if(switchData.LoadInLaserProfile)
                {
                    externalSources.updateLaserPower(timeData.TotalTime);
                }
                plasmaParams.calculatePlasmaFrequency(&variableData.at(lastStepDataIndex), constantData);
                plasmaParams.calculateMinImpactParameter(&variableData.at(lastStepDataIndex), constantData);
                plasmaParams.calculateThermalVelocity(&variableData.at(lastStepDataIndex), constantData);
                plasmaParams.calculateCoulombLogLaser(&variableData.at(lastStepDataIndex), constantData, externalSources.LaserWavelength);
                plasmaParams.calculateCollisionFrequencyEIOverC(&variableData.at(lastStepDataIndex), constantData, true);
                externalSources.calculateBeta(&variableData.at(lastStepDataIndex));
                externalSources.inverseBrem(&variableData.at(lastStepDataIndex));
            }
        }
        if(switchData.ExchangeOn)
        {
            internalSources.exchange(&variableData.at(lastStepDataIndex), constantData, timeData);        
        }
        if(switchData.HeatConductionOn)
        {
            internalSources.calculateKappa(&variableData.at(lastStepDataIndex), constantData); //Spitzer Harm
            if(!switchData.CoupleDivQ)
            {
                internalSources.heatFlowE(&variableData.at(lastStepDataIndex), constantData, switchData);                    
                if(switchData.SNBHeatFlow)
                {
                    variableData.at(lastStepDataIndex).LinearInterpolate();
                    plasmaParams.calculateCoulombLogEI(variableData.at(lastStepDataIndex).InterCoulombLog,
                                            variableData.at(lastStepDataIndex).InterTemperatureE,
                                            variableData.at(lastStepDataIndex).InterNumberDensityE, 
                                            variableData.at(lastStepDataIndex).InterZbar,constantData);
                    // if(constantData.Np > 1)
                    // {
                    //     int finish = allThreads.update(snbSources,  &variableData.at(lastStepDataIndex), 
                    //     constantData);
                    //     if(finish != 0)
                    //     {
                    //     std::cerr << "Possible that all threads have not finished" << std::endl;
                    //     }
                    //     for(auto &i:snbSources)
                    //     {
                        //     i.SNBGather(&variableData.at(lastStepDataIndex));
                        // }
                    //     snbSources.at(0).SNBCorrect(&variableData.at(lastStepDataIndex));
                    // }
                    // else
                    // {
                    for(auto &i:snbSources)
                    {
                        // i.snbHeatFlowCorrection(&variableData.at(lastStepDataIndex), constantData, snbMatrixSolver.at(0));
                        i.snbHeatFlowCorrection(&variableData.at(lastStepDataIndex), constantData);
                        i.SNBGather(&variableData.at(lastStepDataIndex));
                    }
                    snbSources.at(0).SNBCorrect(&variableData.at(lastStepDataIndex));
                    // }
                }
                if(switchData.CoupleOperatorSplit)
                {
                    variableData.at(lastStepDataIndex).CoupleOperatorSplitHeatFlowE = variableData.at(lastStepDataIndex).HeatFlowE; 
                }
                if(switchData.CoupleMulti)
                {
                    variableData.at(lastStepDataIndex).HeatFlowMultiplier = variableData.at(lastStepDataIndex -1 ).HeatFlowMultiplier; // Carry over Multipliers
                    internalSources.multiplierHeatFlowE(&variableData.at(lastStepDataIndex), constantData);                    
                }
                if(switchData.CoupleSubtract)
                {
                    variableData.at(lastStepDataIndex).VFPHeatFlow = variableData.at(lastStepDataIndex -1 ).VFPHeatFlow; // Carry over dq
                    internalSources.subtractHeatFlowE(&variableData.at(lastStepDataIndex));
                }
            }
            if(switchData.CoupleDivQ)
            {
                if(switchData.CoupleOperatorSplit)
                {
                    internalSources.heatFlowE(&variableData.at(lastStepDataIndex), constantData, switchData);                   //Spitzer-Harm 
                    variableData.at(lastStepDataIndex).CoupleOperatorSplitHeatFlowE = variableData.at(lastStepDataIndex).HeatFlowE; //CopeplerOperatorSplitHeatFlow == Spitzer-Harm 
                    internalSources.operatorSplitThermalConducE(&variableData.at(lastStepDataIndex));
                }
                variableData.at(lastStepDataIndex).HeatConductionE = variableData.at(lastStepDataIndex - 1).HeatConductionE; //Load in Div.q 
            }
            else
            {
                if(switchData.CoupleOperatorSplit)
                {
                    internalSources.operatorSplitThermalConducE(&variableData.at(lastStepDataIndex));
                }
                internalSources.thermalConducE(&variableData.at(lastStepDataIndex));
            }
            internalSources.heatFlowI(&variableData.at(lastStepDataIndex), constantData);
            internalSources.thermalConducI(&variableData.at(lastStepDataIndex));
        }

        if(switchData.RadTransportOn)
        {
            int j = 1;
            for(auto &i:radiationTransport)
            {
             
                if(switchData.LeakSource)
                {
                    i.setLeakArea(constantData, timeData.TotalTime);
                }
                i.updateRadiationQuantities(&variableData.at(lastStepDataIndex), constantData, timeData);
                i.gather(&variableData.at(lastStepDataIndex));
                if((constantData.Np == 1)&& (j < constantData.RadNg))
                {
                    radiationTransport.at(j).copySearchEntries(i);
                }
                j++;

            }
            radiationTransport.at(0).updateRadiationTemperature(&variableData.at(lastStepDataIndex), constantData);
            radiationTransport.at(0).getOpacities(&variableData.at(lastStepDataIndex));
        }

        //timeData Step update        
        timeData.setTotalTime();
        timeData.incrementStep();
        if(timeData.Step % constantData.OutputFrequency == 0)
        {
            variableData.back().dump(ioData, timeData.Step);
            timeData.dump(ioData);
        }
        if(switchData.AdapativeTimeStep)
        {
            timeData.calculateNewTime(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex));    
        }
        
        //Printing to terminal for debug. 
        if(verbose_output)
        {
            variableData.back().terminalPrint(); 
        }
        i++;
    }
    //Kill all jobs and join them i.e. kill threads.
    // if((switchData.RadTransportOn)||(switchData.SNBHeatFlow)||(switchData.ImplicitHeatConduction))
    // {
    //     for(auto &i:matrixSolver)
    //     {
    //         i.cleanUpPETScObjects();
        // }
        // PetscFinalize();
    // }
} 
