#include "IntSourceTests.h"

void InitTwoBath(GridData *gridData, FixedData const &fixedData, int nx, double L, double nu,
                double Tu, double Td)
{
    double dx = L / (nx + 1); 
    gridData->TemperatureE.resize(nx, 0);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i * dx);
        gridData->Velocity.push_back(0);
        if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }

    for(int i = 1; i < nx - 1; i++)
    {
        gridData->TemperatureE[i] = pow(pow(Tu, 3.5) - (gridData->CellCenteredCoord[i] / L) 
                                    *(pow(Tu, 3.5) - pow(Td, 3.5)),(1/3.5));
    }
    gridData->TemperatureE[0] = Td;
    gridData->TemperatureE[nx - 1] = Tu;
    gridData->TemperatureI = gridData->TemperatureE;
    for(int i = 0; i < nx; i++)
    {
        gridData->NumberDensityE[i] =  nu * (Tu/gridData->TemperatureE[i]);
        gridData->NumberDensityI[i] = gridData->NumberDensityE[i];
        gridData->Density.push_back(gridData->NumberDensityI[i] * fixedData.PROTON_MASS);
    }
    for(int i = 0; i < nx; i++)
    {
        gridData->Mass.push_back(gridData->Density[i] * (gridData->CellWallCoord[i + 1] - gridData->CellWallCoord[i]));
    }
}

TEST(IntSourceTwoBathTests, DISABLED_ExplicitTwoBathNoChange)
{   
    int nx = 100;
    double L = 1;
    double single_source, Td = 30, Tu = 30000, nu = 1e21;
    TimeStep timeData(nx);
    FixedData fixedData;
    FluidDynamics Hydro(nx);
    IntSource heatflow(nx);
    IdealGasEoS EoS(0, nx, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    fixedData.Ar.resize(nx, 1);
    InitTwoBath(&initGridData, fixedData, nx, L, nu, Tu, Td);
    EoS.updatePressureTerms(&initGridData);
    EoS.updateSpecificHeatCapacity(&initGridData);
    std::fill(initGridData.Zbar.begin(), initGridData.Zbar.end(), 1);
    coulombLog.calculateCoulombLogEI(initGridData.CoulombLogEI, initGridData.TemperatureE,
                                     initGridData.NumberDensityE, initGridData.Zbar, fixedData);
    heatflow.calculateKappa(&initGridData, fixedData);
    heatflow.heatFlowE(&initGridData, fixedData, switchData);
    heatflow.thermalConducE(&initGridData);
    variableData.push_back(initGridData);
    
    timeData.Dt05 = 1E-14;
    timeData.TotalTime = 0;
    timeData.MaxSteps = 1E4;
    fixedData.RightFluidBoundaryCondition = "r";
    switchData.AdibaticModeON = false;
    
    for(int j = 0; j < timeData.MaxSteps; j++)
    {
        
        if(variableData.size() > 4)
        {
                variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        
        //Units being tested
        //Hydro Units
        variableData.at(lastStepDataIndex).Mass =  variableData.at(lastStepDataIndex - 1).Mass; //Mass stays fixed between steps 
        Hydro.updateMomentum(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), fixedData, timeData);
        Hydro.updateCoords(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData);
        Hydro.updateCellCentreCoords(&variableData.at(lastStepDataIndex));
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateDensity(&variableData.at(lastStepDataIndex));
        Hydro.updateNumberDensityI(&variableData.at(lastStepDataIndex), fixedData);
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);
        variableData.at(lastStepDataIndex).TemperatureE[0] = Td;
        variableData.at(lastStepDataIndex).TemperatureE[nx - 1] = Tu;

        variableData.at(lastStepDataIndex).TemperatureI = variableData.at(lastStepDataIndex).TemperatureE;
        
        //Ideal Gas EoS
        std::fill(variableData.at(lastStepDataIndex).Zbar.begin(), variableData.at(lastStepDataIndex).Zbar.end(), 1);
        Hydro.updateNumberDensityE(&variableData.at(lastStepDataIndex));
        EoS.updateSpecificHeatCapacity(&variableData.at(lastStepDataIndex));
        EoS.updatePressureTerms(&variableData.at(lastStepDataIndex));
        coulombLog.calculateCoulombLogEI(variableData.at(lastStepDataIndex).CoulombLogEI,
                                            variableData.at(lastStepDataIndex).TemperatureE,
                                            variableData.at(lastStepDataIndex).NumberDensityE, 
                                            variableData.at(lastStepDataIndex).Zbar,fixedData);
        
        //HeatFlow
        heatflow.calculateKappa(&variableData.at(lastStepDataIndex), fixedData);
        heatflow.heatFlowE(&variableData.at(lastStepDataIndex), fixedData, switchData);
        heatflow.thermalConducE(&variableData.at(lastStepDataIndex));
    }

    std::vector<double> TemperatureChange;
    int lastStepDataIndex = variableData.size();  // index values is -1 of size
    for(int t = 0; t < lastStepDataIndex; t++)
    {
        TemperatureChange.push_back(abs(initGridData.TemperatureE[t] - 
                                variableData.at(lastStepDataIndex - 1).TemperatureE[t]) / initGridData.TemperatureE[t]);
    }
    auto max_error = std::max_element(std::begin(TemperatureChange), std::end(TemperatureChange));
    EXPECT_TRUE(*max_error< 1e-4);
}
TEST(IntSourceTwoBathTests, ImplicitTwoBathNoChange)
{   
    int nx = 100;
    double L = 1;
    double single_source, Td = 30, Tu = 30000, nu = 1e21;
    TimeStep timeData(nx);
    FixedData fixedData;
    FluidDynamics Hydro(nx);
    IntSource heatflow(nx);
    IdealGasEoS EoS(0, nx, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    fixedData.Ar.resize(nx, 1);
    fixedData.heatflxlmE = 1.0;
    InitTwoBath(&initGridData, fixedData, nx, L, nu, Tu, Td);
    EoS.updatePressureTerms(&initGridData);
    EoS.updateSpecificHeatCapacity(&initGridData);
    std::fill(initGridData.Zbar.begin(), initGridData.Zbar.end(), 1);
    coulombLog.calculateCoulombLogEI(initGridData.CoulombLogEI, initGridData.TemperatureE,
                                     initGridData.NumberDensityE, initGridData.Zbar, fixedData);
    heatflow.calculateKappa(&initGridData, fixedData);
    heatflow.heatFlowE(&initGridData, fixedData, switchData);
    heatflow.thermalConducE(&initGridData);
    Hydro.updateDxs(&initGridData);
    variableData.push_back(initGridData);
    
    timeData.Dt05 = 1E-14;
    timeData.TotalTime = 0;
    timeData.MaxSteps = 1E4;
    fixedData.RightFluidBoundaryCondition = "r";
    switchData.AdibaticModeON = false;
    
    for(int j = 0; j < timeData.MaxSteps; j++)
    {
        if(variableData.size() > 4)
        {
                variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        
        //Units being tested
        //Hydro Units
        variableData.at(lastStepDataIndex).Mass =  variableData.at(lastStepDataIndex - 1).Mass; //Mass stays fixed between steps 
        Hydro.updateMomentum(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), fixedData, timeData);
        Hydro.updateCoords(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData);
        Hydro.updateCellCentreCoords(&variableData.at(lastStepDataIndex));
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateDensity(&variableData.at(lastStepDataIndex));
        Hydro.updateNumberDensityI(&variableData.at(lastStepDataIndex), fixedData);
        std::fill(variableData.at(lastStepDataIndex - 1).HeatConductionE.begin(), variableData.at(lastStepDataIndex - 1).HeatConductionE.end(), 0);
        variableData.at(lastStepDataIndex).SpecificHeatE = variableData.at(lastStepDataIndex - 1).SpecificHeatE;
        variableData.at(lastStepDataIndex).HeatKappaE = variableData.at(lastStepDataIndex - 1).HeatKappaE;
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);

        variableData.at(lastStepDataIndex).NumberDensityE = variableData.at(lastStepDataIndex -1 ).NumberDensityE; // Carry over Multipliers
        variableData.at(lastStepDataIndex).Zbar = variableData.at(lastStepDataIndex -1 ).Zbar; // Carry over Multipliers
        coulombLog.calculateCoulombLogEI(variableData.at(lastStepDataIndex).CoulombLogEI,
                                        variableData.at(lastStepDataIndex).TemperatureE,
                                        variableData.at(lastStepDataIndex).NumberDensityE, 
                                        variableData.at(lastStepDataIndex).Zbar,fixedData);
        heatflow.calculateKappa(&variableData.at(lastStepDataIndex), fixedData);
        heatflow.calculateImplicitKappaE(&variableData.at(lastStepDataIndex), fixedData);
        // Hydro.updateTemperatureEImplicitHeatConduction(&variableData.at(lastStepDataIndex), timeData, switchData);
        Hydro.updateTemperatureImplicitHeatConduction(variableData.at(lastStepDataIndex).TemperatureE, variableData.at(lastStepDataIndex).HeatKappaE, variableData.at(lastStepDataIndex).SpecificHeatE,
                                            variableData.at(lastStepDataIndex).Mass, variableData.at(lastStepDataIndex).Dx_k, timeData.Dt05);
        variableData.at(lastStepDataIndex).TemperatureE[0] = Td;
        variableData.at(lastStepDataIndex).TemperatureE[nx - 1] = Tu;

        variableData.at(lastStepDataIndex).TemperatureI = variableData.at(lastStepDataIndex).TemperatureE;
        
        //Ideal Gas EoS
        std::fill(variableData.at(lastStepDataIndex).Zbar.begin(), variableData.at(lastStepDataIndex).Zbar.end(), 1);
        Hydro.updateNumberDensityE(&variableData.at(lastStepDataIndex));
        EoS.updateSpecificHeatCapacity(&variableData.at(lastStepDataIndex));
        EoS.updatePressureTerms(&variableData.at(lastStepDataIndex));
        coulombLog.calculateCoulombLogEI(variableData.at(lastStepDataIndex).CoulombLogEI,
                                            variableData.at(lastStepDataIndex).TemperatureE,
                                            variableData.at(lastStepDataIndex).NumberDensityE, 
                                            variableData.at(lastStepDataIndex).Zbar,fixedData);
        
        //HeatFlow
        heatflow.calculateKappa(&variableData.at(lastStepDataIndex), fixedData);
        heatflow.heatFlowE(&variableData.at(lastStepDataIndex), fixedData, switchData);
        heatflow.thermalConducE(&variableData.at(lastStepDataIndex));
    }

    std::vector<double> TemperatureChange;
    int lastStepDataIndex = variableData.size();  // index values is -1 of size
    for(int t = 0; t < lastStepDataIndex; t++)
    {
        TemperatureChange.push_back(abs(initGridData.TemperatureE[t] - 
                                variableData.at(lastStepDataIndex - 1).TemperatureE[t]) / initGridData.TemperatureE[t]);
    }
    auto max_error = std::max_element(std::begin(TemperatureChange), std::end(TemperatureChange));
    EXPECT_TRUE(*max_error< 1e-4);
}
TEST(IntSourceNonLocalHeatRecreationTest, subtractCheck)
{
        
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/VFP_SUBTRACT/";
    int nx = 200;
    FixedData fixedData;
    fixedData.initVectors(localPath, true);
    IntSource heatflow(nx);
    FluidDynamics Hydro(nx);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    switchData.Couple = true;
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    initGridData.init(localPath, false, true, 3);
    Hydro.updateNumberDensityI(&initGridData, fixedData);
    Hydro.updateCellCentreCoords(&initGridData);
    Hydro.updateDxs(&initGridData);
    initGridData.Zbar = fixedData.Z;
    Hydro.updateNumberDensityE(&initGridData);
    // coulombLog.calculateCoulombLogEI(initGridData.CoulombLogEI, initGridData.TemperatureE,
    //                                  initGridData.NumberDensityE, initGridData.Zbar, fixedData);
    coulombLog.calculateCoulombLogEI(&initGridData, fixedData);
    heatflow.calculateKappa(&initGridData, fixedData);
    heatflow.heatFlowE(&initGridData, fixedData, switchData);
    heatflow.subtractHeatFlowE(&initGridData);
    std::vector<double> relative_error;
    for(int i = 0; i <= nx; i++)
    {
        double error = abs(abs(initGridData.HeatFlowE[i]) - abs(True_heat_flow_value[i]))
                        /abs(True_heat_flow_value[i]);
        if(std::isnan(error))
        {
            continue;
        }
        relative_error.push_back(error);
    }
    auto max_error = std::max_element(std::begin(relative_error), std::end(relative_error));
    EXPECT_TRUE(*max_error< 1e-2);
}
TEST(IntSourceNonLocalHeatRecreationTest, divQCheck)
{
        
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/VFP_DivQ/";
    int nx = 200;
    FixedData fixedData;
    fixedData.initVectors(localPath, true);
    IntSource heatflow(nx);
    FluidDynamics Hydro(nx);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    initGridData.init(localPath, false, true, 2);
    Hydro.updateNumberDensityI(&initGridData, fixedData);
    Hydro.updateCellCentreCoords(&initGridData);
    Hydro.updateDxs(&initGridData);
    initGridData.Zbar = fixedData.Z;
    Hydro.updateNumberDensityE(&initGridData);
    std::vector<double> relative_error;
    for(int i = 0; i <= nx; i++)
    {
        double error;
        if (abs(initGridData.LoadInDivQ[i]) < 1e-30)
        {
            error = 1e-10;
        }
        else
        {
            error = abs(abs(true_load_div_q[i]) - abs(initGridData.LoadInDivQ[i]))
                            /abs(true_load_div_q[i]);
        }
        if(std::isnan(error))
        {
            continue;
        }
        relative_error.push_back(error);
    }
    auto max_error = std::max_element(std::begin(relative_error), std::end(relative_error));
    EXPECT_TRUE(*max_error< 1e-9);
}
// TEST(IntSourceNonLocalHeatRecreationTest, SNBMultiCheck)
// {
        
//     std::string loadPath = getenv("TEST_PATH");
//     std::string localPath = loadPath + "/SNB_MULTIPLIER_RECREATER_CHECK/";
//     int nx = 200;
//     FixedData fixedData;
//     fixedData.initVectors(localPath, true);
//     IntSource heatflow(nx);
//     FluidDynamics Hydro(nx);
//     PlasmaParameters coulombLog(nx);
//     Switches switchData;
//     switchData.SNBHeatFlow = true;
//     std::vector<GridData> variableData;
//     GridData initGridData(nx, true);
//     MatrixSolver matrixSolver(nx);
//     matrixSolver.createPETScObjects();
//     initGridData.init(localPath, false, true, 1);
//     Hydro.updateNumberDensityI(&initGridData, fixedData);
//     Hydro.updateCellCentreCoords(&initGridData);
//     Hydro.updateDxs(&initGridData);
//     initGridData.Zbar = fixedData.Z;
//     Hydro.updateNumberDensityE(&initGridData);
//     coulombLog.calculateCoulombLogEI(&initGridData, fixedData);
//     auto max_T = *std::max_element(std::begin(initGridData.TemperatureE), std::end(initGridData.TemperatureE));
//     heatflow.InitSNB(20, 12*max_T);
//     heatflow.calculateKappa(&initGridData, fixedData);
//     heatflow.heatFlowE(&initGridData, fixedData, switchData);
//     initGridData.LinearInterpolate();
//     coulombLog.calculateCoulombLogEI(initGridData.InterCoulombLog,
//                             initGridData.InterTemperatureE,
//                             initGridData.InterNumberDensityE, 
//                             initGridData.InterZbar,fixedData);
//     // internalSources.createPETScObjects();
//     heatflow.snbHeatFlowCorrection(&initGridData, fixedData, matrixSolver);
//     std::vector<double> relative_error;
//     double error;
//     for(int i = 0; i <= nx; i++)
//     {
//         error = abs(abs(True_heat_flow_value[i]) - abs(initGridData.HeatFlowE[i]))
//                         /abs(True_heat_flow_value[i]);
//         if(std::isnan(error))
//         {
//             continue;
//         }
//         relative_error.push_back(error);
//     }
//     auto max_error = std::max_element(std::begin(relative_error), std::end(relative_error));
//     EXPECT_TRUE(*max_error< 1e-4);
// }
// TEST(IntSourceNonLocalHeatRecreationTest, DISABLED_SNBSubtractCheck)
// {
        
//     std::string loadPath = getenv("TEST_PATH");
//     std::string localPath = loadPath + "/SNB_SUBTRACT_RECREATER_CHECK/";
//     int nx = 200;
//     FixedData fixedData;
//     fixedData.initVectors(localPath, true);
//     IntSource heatflow(nx);
//     FluidDynamics Hydro(nx);
//     PlasmaParameters coulombLog(nx);
//     Switches switchData;
//     switchData.SNBHeatFlow = true;
//     std::vector<GridData> variableData;
//     GridData initGridData(nx, true);
//     MatrixSolver matrixSolver(nx);
//     matrixSolver.createPETScObjects();
//     initGridData.init(localPath, false, true, 3);
//     Hydro.updateNumberDensityI(&initGridData, fixedData);
//     Hydro.updateCellCentreCoords(&initGridData);
//     Hydro.updateDxs(&initGridData);
//     initGridData.Zbar = fixedData.Z;
//     Hydro.updateNumberDensityE(&initGridData);
//     coulombLog.calculateCoulombLogEI(&initGridData, fixedData);
//     auto max_T = *std::max_element(std::begin(initGridData.TemperatureE), std::end(initGridData.TemperatureE));
//     heatflow.InitSNB(20, 12*max_T);
//     heatflow.calculateKappa(&initGridData, fixedData);
//     heatflow.heatFlowE(&initGridData, fixedData, switchData);
//     initGridData.LinearInterpolate();
//     coulombLog.calculateCoulombLogEI(initGridData.InterCoulombLog,
//                             initGridData.InterTemperatureE,
//                             initGridData.InterNumberDensityE, 
//                             initGridData.InterZbar,fixedData);
//     // internalSources.createPETScObjects();
//     heatflow.snbHeatFlowCorrection(&initGridData, fixedData, matrixSolver);
//     std::vector<double> relative_error;
//     double error;
//     for(int i = 0; i <= nx; i++)
//     {
//         error = abs(abs(True_heat_flow_value[i]) - abs(initGridData.HeatFlowE[i]))
//                         /abs(True_heat_flow_value[i]);
//         if(std::isnan(error))
//         {
//             continue;
//         }
//         relative_error.push_back(error);
//     }
//     auto max_error = std::max_element(std::begin(relative_error), std::end(relative_error));
//     EXPECT_TRUE(*max_error< 1e-4);
// }

void loadInBrodrickQuants(std::string localPath, GridData &gridData, std::vector<double> &true_qe, bool seperated)
{
    std::vector<std::string> vars;
    if(seperated)
    {
        vars = {"coord", "Te", "Z", "ne",
                    "qe_sep"};
    }
    else
    {
        vars = {"coord", "Te", "Z", "ne",
                    "qe_avg"};
    }
    int no_var = vars.size();
    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::string file_path = localPath + "/" +vars[i] + ".txt";
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;
        if (initFile.is_open())
        {
            int j = 0;
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    gridData.CellWallCoord[j] = stod(strInput);
                }
                else if( i==1)
                {
                    gridData.TemperatureE[j] = stod(strInput);

                }
                else if( i==2)
                {

                    gridData.Zbar[j] = stod(strInput);
                }
                else if( i==3)
                {
                    gridData.NumberDensityE[j] = stod(strInput);

                }
                else if( i==4)
                {   
                    true_qe[j] = stod(strInput) * -1e4;
                }
                j++;
            }
        }
    }
}
TEST(IntSourceSNBCheck, SNBBrodrickSeperatedRecreate)
{
        
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/SNB_BRODRICK_RECREATE/";
    int nx = 600;
    FixedData fixedData;
    fixedData.Ng = 50;
    fixedData.Nx = nx;
    IntSource heatflow(nx);
    FluidDynamics Hydro(nx);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    switchData.SNBHeatFlow = true;
    std::vector<GridData> variableData;
    GridData initGridData(nx, false);
    std::vector<SNBSource> snbSources(fixedData.Ng, fixedData.Nx);
    std::vector<double> true_qe(nx + 1,0.0);
    loadInBrodrickQuants(localPath, initGridData, true_qe, true);
    Hydro.updateCellCentreCoords(&initGridData);
    Hydro.updateDxs(&initGridData);
    // coulombLog.calculateCoulombLogEI(&initGridData, fixedData);
    std::vector<double> Ar(nx, 10);
    int j = 0;
    for(auto &i:snbSources)
    {
        i.InitSNB(fixedData.Ng, j,
                *max_element(initGridData.TemperatureE.begin(),
                initGridData.TemperatureE.end())*20.0);
        j++;
    }
    std::fill(initGridData.CoulombLogEI.begin(), initGridData.CoulombLogEI.end(), 2.1484);
    std::fill(initGridData.TemperatureI.begin(), initGridData.TemperatureI.end(), 2.1484);
    fixedData.Ar = Ar;
    heatflow.calculateKappa(&initGridData, fixedData);
    heatflow.heatFlowE(&initGridData, fixedData, switchData);
    initGridData.LinearInterpolate();
    std::fill(initGridData.InterCoulombLog.begin(), initGridData.InterCoulombLog.end(), 2.1484);

    for(auto &i:snbSources)
    {
        i.snbHeatFlowCorrection(&initGridData, fixedData);
        i.SNBGather(&initGridData);
        i.SNBGatherH(&initGridData);
    }
    snbSources.at(0).SNBCorrect(&initGridData);
    // initGridData.dump(ioData, 0);
    std::vector<double> relative_error;
    for(int i = 0; i <= nx; i++)
    {
        double error = abs(abs(initGridData.HeatFlowE[i]) - abs(true_qe[i]))
                        /abs(true_qe[i]);
        if(std::isnan(error))
        {
            relative_error.push_back(0.0);
            continue;
        }
        relative_error.push_back(error);
    }
    double max_error = *std::max_element(std::begin(relative_error), std::end(relative_error));
    int arg_index = std::distance(relative_error.begin(),std::max_element(relative_error.begin(), relative_error.end()));
    EXPECT_TRUE(max_error< 1e-2);
}

TEST(IntSourceSNBCheck, SNBBrodrickAveragedRecreate)
{
        
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/SNB_BRODRICK_RECREATE/";
    
    int nx = 600;
    FixedData fixedData;
    fixedData.Ng = 50;
    fixedData.Nx = nx;
    IntSource heatflow(nx);
    FluidDynamics Hydro(nx);
    PlasmaParameters coulombLog(nx);
    Switches switchData;
    switchData.SNBHeatFlow = true;
    std::vector<GridData> variableData;
    GridData initGridData(nx, false);
    std::vector<SNBSource> snbSources(fixedData.Ng, fixedData.Nx);
    std::vector<double> true_qe(nx + 1,0.0);
    loadInBrodrickQuants(localPath, initGridData, true_qe, false);
    Hydro.updateCellCentreCoords(&initGridData);
    Hydro.updateDxs(&initGridData);
    // coulombLog.calculateCoulombLogEI(&initGridData, fixedData);
    std::vector<double> Ar(nx, 10);
    int j = 0;
    for(auto &i:snbSources)
    {
        i.InitSNB(fixedData.Ng, j,
                *max_element(initGridData.TemperatureE.begin(),
                initGridData.TemperatureE.end())*20.0);
        i.switchToAveragedSNB();
        j++;
    }
    std::fill(initGridData.CoulombLogEI.begin(), initGridData.CoulombLogEI.end(), 2.1484);
    std::fill(initGridData.TemperatureI.begin(), initGridData.TemperatureI.end(), 2.1484);
    fixedData.Ar = Ar;
    heatflow.calculateKappa(&initGridData, fixedData);
    heatflow.heatFlowE(&initGridData, fixedData, switchData);
    initGridData.LinearInterpolate();
    std::fill(initGridData.InterCoulombLog.begin(), initGridData.InterCoulombLog.end(), 2.1484);

    for(auto &i:snbSources)
    {
        i.snbHeatFlowCorrection(&initGridData, fixedData);
        i.SNBGather(&initGridData);
        i.SNBGatherH(&initGridData);
    }
    snbSources.at(0).SNBCorrect(&initGridData);
    // initGridData.dump(ioData, 0);
    std::vector<double> relative_error;
    for(int i = 0; i <= nx; i++)
    {
        double error = abs(abs(initGridData.HeatFlowE[i]) - abs(true_qe[i]))
                        /abs(true_qe[i]);
        if(std::isnan(error))
        {
            relative_error.push_back(0.0);
            continue;
        }
        relative_error.push_back(error);
    }
    double max_error = *std::max_element(std::begin(relative_error), std::end(relative_error));
    int arg_index = std::distance(relative_error.begin(),std::max_element(relative_error.begin(), relative_error.end()));
    EXPECT_TRUE(max_error< 1e-2);
}

TEST(IntSourceHeatExchange, HeatExchangeTestConstTi)
{
    int nx = 1;
    FixedData fixedData;
    fixedData.Nx = nx;
    IntSource exchange(nx);
    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    TimeStep timeData(nx);
    timeData.MaxSteps = 10000;
    timeData.Dt05 = 1e-9;
    timeData.Dt1 = 1e-9;
    double Te_0 = 100 * fixedData.ev_to_k;
    double Ti_0 = 100000 * fixedData.ev_to_k;
    double specificHeatCapacity = 50;
    initGridData.TemperatureE.resize(nx, Te_0);
    initGridData.TemperatureI.resize(nx, Ti_0);
    initGridData.SpecificHeatE[0] = specificHeatCapacity;
    initGridData.Density.resize(nx, 1);
    initGridData.DpDtE[0] = 1;
    double ex_coef = 4e+5; 
    exchange.setExchangeCoefficient(&initGridData, ex_coef);
    variableData.push_back(initGridData);

    
    auto expected_result = [](double T0, double Ti, double Cv, double ex_coef, double t){return (Ti - (Ti - T0) *exp(-1*(ex_coef * t) / Cv) );};

    for(int j = 0; j < timeData.MaxSteps; j++)
    {
        if(variableData.size() > 4)
        {
            variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        
        //Units being tested
        //Hydro Units
        
        variableData.at(lastStepDataIndex).Density = variableData.at(lastStepDataIndex - 1).Density;
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);
        variableData.at(lastStepDataIndex).SpecificHeatE = variableData.at(lastStepDataIndex - 1).SpecificHeatE;
        variableData.at(lastStepDataIndex).TemperatureI = variableData.at(lastStepDataIndex - 1).TemperatureI;
        variableData.at(lastStepDataIndex).DpDtE = variableData.at(lastStepDataIndex - 1).DpDtE;
        exchange.setExchangeCoefficient(&variableData.at(lastStepDataIndex), ex_coef);
        timeData.incrementStep();
        timeData.setTotalTime();
        double analytical_Te = expected_result(Te_0, Ti_0, specificHeatCapacity, ex_coef, timeData.TotalTime);
        double relative_error = ((variableData.at(lastStepDataIndex).TemperatureE[0] - analytical_Te) / analytical_Te);
        EXPECT_TRUE(relative_error < 1e-2);
    }
}

TEST(IntSourceHeatExchange, HeatExchangeTestConstTe)
{

    int nx = 1;
    FixedData fixedData;
    fixedData.Nx = nx;
    IntSource exchange(nx);
    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    TimeStep timeData(nx);
    PlasmaParameters plasmaParams(nx);
    timeData.MaxSteps = 10000;
    timeData.Dt05 = 1e-15;
    timeData.Dt1 = 1e-15;
    double Te_0 = 100 * fixedData.ev_to_k;
    double Ti_0 = 100000 * fixedData.ev_to_k;
    double specificHeatCapacity = 50;
    initGridData.TemperatureE.resize(nx, Te_0);
    initGridData.TemperatureI.resize(nx, Ti_0);
    initGridData.SpecificHeatI[0] = specificHeatCapacity;
    initGridData.SpecificHeatE[0] = specificHeatCapacity;
    initGridData.Density.resize(nx, 1);
    initGridData.DpDtI[0] = 1;
    initGridData.Zbar[0] = 2;
    fixedData.Ar.resize(nx, 4);
    initGridData.CoulombLogEI[0] = 10;
    initGridData.NumberDensityI[0] = initGridData.Density[0] / (fixedData.Ar[0] * fixedData.PROTON_MASS);
    initGridData.NumberDensityE[0] = initGridData.NumberDensityI[0] * initGridData.Zbar[0];
    plasmaParams.calculateCollisionFrequencyEI(&initGridData, fixedData, false);
    exchange.exchange(&initGridData, fixedData, timeData);
    variableData.push_back(initGridData);
    auto expected_result = [](double T0, double Ti, double Cv, double ex_coef, double t){return (T0 + (Ti - T0) *exp(-1*(ex_coef * t) / Cv) );};

    for(int j = 0; j < timeData.MaxSteps; j++)
    {
        if(variableData.size() > 4)
        {
            variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        
        //Units being tested
        //Hydro Units
        
        variableData.at(lastStepDataIndex).Density = variableData.at(lastStepDataIndex - 1).Density;
        variableData.at(lastStepDataIndex).NumberDensityI = variableData.at(lastStepDataIndex - 1).NumberDensityI;
        Hydro.updateTemperatureI(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);
        variableData.at(lastStepDataIndex).SpecificHeatI = variableData.at(lastStepDataIndex - 1).SpecificHeatI;
        variableData.at(lastStepDataIndex).SpecificHeatE = variableData.at(lastStepDataIndex - 1).SpecificHeatE;
        variableData.at(lastStepDataIndex).NumberDensityE = variableData.at(lastStepDataIndex - 1).NumberDensityE;
        variableData.at(lastStepDataIndex).TemperatureE = variableData.at(lastStepDataIndex - 1).TemperatureE;
        variableData.at(lastStepDataIndex).DpDtI = variableData.at(lastStepDataIndex - 1).DpDtI;
        variableData.at(lastStepDataIndex).CoulombLogEI = variableData.at(lastStepDataIndex - 1).CoulombLogEI;
        variableData.at(lastStepDataIndex).Zbar = variableData.at(lastStepDataIndex - 1).Zbar;
        double standard_Exchange = 4.9044757184917434e-05 * variableData.at(lastStepDataIndex).Zbar[0] * 
                                    variableData.at(lastStepDataIndex).Zbar[0] * (1.0 /(fixedData.Ar[0] * fixedData.Ar[0])) * 
                                    variableData.at(lastStepDataIndex).CoulombLogEI[0] * variableData.at(lastStepDataIndex).NumberDensityE[0] * 
                                    pow(variableData.at(lastStepDataIndex).TemperatureE[0], -1.5);         
        plasmaParams.calculateCollisionFrequencyEI(&variableData.at(lastStepDataIndex), fixedData, false);
        exchange.exchange(&variableData.at(lastStepDataIndex), fixedData, timeData);
        timeData.incrementStep();
        timeData.setTotalTime();
        double analytical_Ti = expected_result(Te_0, Ti_0, specificHeatCapacity, standard_Exchange, timeData.TotalTime);
        double relative_error = ((variableData.at(lastStepDataIndex).TemperatureI[0] - analytical_Ti) / analytical_Ti);
        EXPECT_TRUE(abs(relative_error) < 1e-2);
    }
}


// INSTANTIATE_TEST_SUITE_P(
//     TwoBathProblem,
//     IntSourceTwoBathTests,
//     ::testing::Values(
//         std::make_tuple(100, 1, 1E-17),
        // std::make_tuple(100, 1, 1E-15)));
