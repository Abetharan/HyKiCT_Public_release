#include "ImplicitRadiationTests.h"
void InitImplicitSuOlsen(GridData *gridData, FixedData const &fixedData, int nx)
{
    
    double dx = 10E-2 / (nx + 1);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i * dx);
        // if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }
    for(int i = 0; i < nx; i++)
    {
        gridData->Zbar[i] = 1;
        gridData->RadEnergyDensity[i] = 0;
        gridData->RadRossAbsorptionOpacity[i] = 1E-1;
        gridData->RadPlanckAbsorptionOpacity[i] = 1E-1;
        gridData->RadPlanckEmissionOpacity[i] = 1E-1;
        gridData->Density.push_back(1000);
        gridData->Mass.push_back(gridData->Density[i] * (gridData->CellWallCoord[i + 1] - gridData->CellWallCoord[i]));
        gridData->TemperatureE.push_back(5e2);
        gridData->TemperatureI.push_back(5e2);
        gridData->SpecificHeatE[i] = (40/gridData->Density[i]) * fixedData.radiationConstant * pow(gridData->TemperatureE[i], 3);
    }
}
void InitImplicitMultiSuOlsen(GridData *gridData, FixedData const &fixedData, int nx, double kappa_bar,double alpha, double density)
{

    double x_max = 102.4 / kappa_bar;
    double dx = x_max / (nx + 1);
    gridData->CellWallCoord.resize(nx+1);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord[i] = i * dx;
        if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }


    for(int i = 0; i < nx; i++)
    {
        // gridData->CellCenteredCoord[i] = (gridData->CellWallCoord[i + 1] + gridData->CellWallCoord[i]) / 2;
        gridData->Zbar[i] = 1;
        gridData->RadEnergyDensity[i] = 0;
        gridData->Density.push_back(density);
        gridData->Mass.push_back(gridData->Density[i] * (gridData->CellWallCoord[i + 1] - gridData->CellWallCoord[i]));
        gridData->TemperatureE.push_back(5e2);
        gridData->TemperatureI.push_back(5e2);
        gridData->SpecificHeatE[i] = (alpha/gridData->Density[i]) * fixedData.radiationConstant * pow(gridData->TemperatureE[i], 3);
    }
}
std::vector<double> ImplicitMultiSuOlsenPlanckIntegral(GridData *gridData, FixedData const& fixedData)
{
    std::vector<double> mPlanck;
    mPlanck.resize(fixedData.Nx);
    for(int i = 0; i < fixedData.Nx; i++)
    {
        mPlanck[i] = (0.5 * (fixedData.radiationConstant * fixedData.SPEED_OF_LIGHT)/(4*M_PI))*pow(gridData->TemperatureE[i], 4);
    }
    return(mPlanck);
}
void ImplicitEquilibrium(GridData *gridData, FixedData const &fixedData, double E_0)
{
    gridData->CellWallCoord.push_back(0);
    gridData->CellCenteredCoord[0] = 0;
    gridData->Density.push_back(25);
    gridData->TemperatureE.push_back(1.1E6);
    gridData->Zbar.push_back(1);
    gridData->SpecificHeatE[0] = 5;
}
void ImplicitTransport(GridData *gridData, FixedData const &fixedData, TimeStep &timeStep,  int nx)
{

    double dx = 1E-2 / (nx + 1);
    timeStep.Dt05 = dx/ (10 * fixedData.SPEED_OF_LIGHT);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i * dx);
        if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }

    for(int i = 0; i < nx; i++)
    {
        gridData->RadEnergyDensity[i] = 1E9;
        gridData->RadRossAbsorptionOpacity[i] = 0.004;
        gridData->RadPlanckAbsorptionOpacity[i] = 0.004;
        gridData->RadPlanckEmissionOpacity[i] = 0.004;
        gridData->Density.push_back(250);
        gridData->TemperatureE.push_back(1.1E6);
    }
}

void ImplicitInitRadSources(GridData *gridData, double T_0, double R, int nx)
{
    double dx = R / (nx + 1);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i * dx);
        if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }

    for(int i = 0; i < nx; i++)
    {
        gridData->TemperatureE.push_back(T_0 * (1 + pow(sin(gridData->CellCenteredCoord[i] / R), 2)));
        gridData->Density.push_back(1.67);
        gridData->SpecificHeatE[i] = 500;
    }
}
void ImplicitInitRadDiffusion(GridData *gridData, double U, double D_0, double R, int nx)
{

    double dx = R / (nx + 1);
    for(int i = 0; i < nx + 1; i++)
    {
        gridData->CellWallCoord.push_back(i * dx);
        if(i > 0)
        {
            gridData->CellCenteredCoord[i - 1] = (gridData->CellWallCoord[i] + gridData->CellWallCoord[i - 1]) / 2;
        }
    }

    for(int i = 0; i < nx; i++)
    {
        if((i !=0) || (i != nx - 1))
        {
            gridData->RadEnergyDensity[i] = U * (1 + pow(sin(gridData->CellCenteredCoord[i] / R), 2));
        }
        gridData->RadDiffusionCoefficient[i] = D_0 * (1 + gridData->CellWallCoord[i] / R); 
    }
}

TEST(ImplicitRadiationTransport, HandlesSources)
{ 
    int nx = 100;
    double single_source, T_0 = 4.81, R = 1E-2, Kappa =4E-3, U =1E10;
    TimeStep timeData(nx);
    timeData.Dt05 = 1E-15;
    timeData.TotalTime = 0;
    timeData.MaxSteps = 100;
    FixedData fixedData;

    fixedData.LeftRadBoundaryCondition = "c";
    fixedData.RightRadBoundaryCondition = "c";
    MultiRadTrans radiationTransport(nx);
    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    ImplicitInitRadSources(&initGridData, T_0, R, nx);
    Hydro.updateDxs(&initGridData);
    radiationTransport.setPlanckianOpacity(Kappa);
    radiationTransport.setRadEnergyDensity(U);
    std::vector<LoadInTables> opacityTable = {LoadInTables{0, fixedData.Nx}};
    opacityTable.at(0).PhotonGrid = {1.00000000e-050, 1.0e+050};
    int j = 0;
    radiationTransport.storeOpacityTables(opacityTable);
    radiationTransport.setEnergyGroup(j);
    radiationTransport.planckIntegral(&initGridData, fixedData);
    radiationTransport.calcFreeFreeEmission(&initGridData, fixedData);
    radiationTransport.fullyImplicit(&initGridData, fixedData, timeData);
    radiationTransport.calcFreeFreeAbsorption(&initGridData, fixedData);
    // radiationTransport.updateRadiationEnergyDensity(&initGridData, fixedData, timeData);
    radiationTransport.gather(&initGridData);
    variableData.push_back(initGridData);

    for(int j = 0; j < timeData.MaxSteps; j++)
    {
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        variableData[lastStepDataIndex] = variableData[lastStepDataIndex - 1];
        //Units being tested
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex- 1), &variableData.at(lastStepDataIndex), timeData, false);
        std::fill(variableData[lastStepDataIndex].RadFFEmission.begin(), variableData[lastStepDataIndex].RadFFEmission.end(), 0.0);
        std::fill(variableData[lastStepDataIndex].RadFFAbsorb.begin(), variableData[lastStepDataIndex].RadFFAbsorb.end(), 0.0);
        std::fill(variableData[lastStepDataIndex].RadEnergyDensity.begin(), variableData[lastStepDataIndex].RadEnergyDensity.end(), 0.0);
        radiationTransport.planckIntegral(&variableData[lastStepDataIndex], fixedData);
        radiationTransport.calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData);
        radiationTransport.calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        // radiationTransport.updateRadiationEnergyDensity(&variableData[lastStepDataIndex], fixedData, timeData);
        radiationTransport.gather(&variableData.at(lastStepDataIndex));
        //Testing 
        std::vector<double> source;
        double expected_value;
        std::vector<double> relative_errors;
        GridData *gridData = &variableData.at(lastStepDataIndex); 
        for(unsigned int i = 0; i < gridData->TemperatureE.size(); i++)
        {
            single_source = -1 * gridData->SpecificHeatE[i] * T_0 * (1 + pow(sin(gridData->CellCenteredCoord[i] / R), 2)) * 
                                    exp(-1*timeData.TotalTime)  - fixedData.SPEED_OF_LIGHT * Kappa * U - 4 * fixedData.STEFAN_BOLTZMANN_CONSTANT * Kappa * 
                                    pow((T_0 * exp(-1 * timeData.TotalTime) * (1 + pow(sin(gridData->CellCenteredCoord[i]/R), 2))), 4);
            gridData->TemperatureE[i] += (1/gridData->SpecificHeatE[i]) *timeData.Dt05 * single_source;
            expected_value = T_0 * (1 + pow(sin(gridData->CellCenteredCoord[i] / R), 2) * exp(-1 * timeData.TotalTime));
            relative_errors.push_back(abs(gridData->TemperatureE[i] - expected_value) / expected_value);
        }
        auto max_error = std::max_element(std::begin(relative_errors), std::end(relative_errors));
        radiationTransport.setRadEnergyDensity(U);
        timeData.setTotalTime();
        ASSERT_LT(*max_error,  1e-8);
    }
}

TEST(ImplicitRadiationTransport, HandlesDiffusion)
{
    //Initialise PETSC 
    // PetscInitializeNoArguments();
    //Init variables
    int nx = 100;
    double single_source, t_0 = 1e-8, R = 1E-2, D_0 =4E-3, U =1E10;
    TimeStep timeData(nx);
    timeData.Dt05 = 1E-17;
    timeData.TotalTime = 0;
    timeData.Tmax = 1E-15;
    FixedData fixedData;
    fixedData.RightRadBoundaryCondition = "d";
    fixedData.LeftRadBoundaryCondition = "d";
    MultiRadTrans radiationTransport(nx);
    FluidDynamics Hydro(nx);
    // MatrixSolver matrixSolver(nx);
    // matrixSolver.createPETScObjects();
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    ImplicitInitRadDiffusion(&initGridData, U, D_0, R, nx);
    initGridData.RadEnergyDensity[0] = U *(1 + pow(sin(initGridData.CellCenteredCoord[0] /R), 2))
                                                *exp((timeData.TotalTime * (t_0 - timeData.TotalTime))/ t_0);

    initGridData.RadEnergyDensity[nx - 1]= U *(1 + pow(sin(initGridData.CellCenteredCoord[nx - 1] /R), 2))
                                                *exp((timeData.TotalTime * (t_0 - timeData.TotalTime))/ t_0);
    radiationTransport.setRadEnergyDensity(initGridData.RadEnergyDensity); 
    radiationTransport.setDiffusionCoefficient(initGridData.RadDiffusionCoefficient);
    fixedData.LeftDirichletEnergyDensity = initGridData.RadEnergyDensity[0];
    fixedData.RightDirichletEnergyDensity = initGridData.RadEnergyDensity[nx - 1];
    Hydro.updateDxs(&initGridData);
    variableData.push_back(initGridData);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {
            break;
        }
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        variableData[lastStepDataIndex] = variableData[lastStepDataIndex - 1];
        if(timeData.TotalTime > 0.)
        {
            fixedData.LeftDirichletEnergyDensity = U *(1 + pow(sin(initGridData.CellCenteredCoord[0] /R), 2))
                                                *exp((timeData.TotalTime * (t_0 - timeData.TotalTime))/ t_0);

            fixedData.RightDirichletEnergyDensity  = U *(1 + pow(sin(initGridData.CellCenteredCoord[nx - 1] /R), 2))
                                                *exp((timeData.TotalTime * (t_0 - timeData.TotalTime))/ t_0);
        }
        GridData *gridData = &variableData.at(lastStepDataIndex); 
        std::fill(gridData->RadEnergyDensity.begin(), gridData->RadEnergyDensity.end(), 0); 
        //Units being tested
        // radiationTransport.ficksDiffusion(gridData, fixedData, timeData);
        radiationTransport.fullyImplicit(gridData, fixedData, timeData);
        // radiationTransport.updateRadiationEnergyDensity(gridData);
        radiationTransport.gather(gridData);

        //Testing
        double expected_value;
        std::vector<double> relative_errors; 
        for(unsigned int i = 0; i < nx; i++)
        {
            double trig = (
                        (1 - 2*(timeData.TotalTime / t_0)) * (1 + pow(sin(gridData->CellCenteredCoord[i] / R), 2))
                        - ((2 *D_0)/pow(R,2)) * ((1 + gridData->CellCenteredCoord[i] / R) * cos((2 * gridData->CellCenteredCoord[i])/R))
                        + sin(gridData->CellCenteredCoord[i]/R) * cos(gridData->CellCenteredCoord[i]/R)
                        );
            single_source = U * exp((timeData.TotalTime * (t_0 - timeData.TotalTime))/ t_0) * trig;  
            gridData->RadEnergyDensity[i] += timeData.Dt05 * single_source;
            expected_value = U * (1 + pow(sin(gridData->CellCenteredCoord[i]/ R), 2)) * 
                                exp((timeData.TotalTime*(t_0 - timeData.TotalTime))/ t_0);
            relative_errors.push_back(abs(gridData->RadEnergyDensity[i] - expected_value) / expected_value);
        }
        timeData.setTotalTime();
        auto max_error = std::max_element(std::begin(relative_errors), std::end(relative_errors));
        ASSERT_LT(*max_error, 1e-8);
    }
}

std::vector<double> ImplicitAnalyticEquilibriumResult(double E_0, double dt, int t_step)
{
    std::vector<double> E_r;
    double T = 1.1E6, rho = 25, kappa_p = 0.04;
    double sigma_b = 5.670374419E-8;
    double c = 3E8;
    for(int i = 0; i <t_step; i++)
    {
        double t = i * dt;
        double B = ((4 * sigma_b) / c) * pow(T, 4);
        E_r.push_back(B - (B - E_0) * exp((-rho*kappa_p*c*t)));
    }
    return(E_r);
}
TEST_P(ImplicitGrayRadiationTransportEquil, Equilibrium_Test)
{

    int nx = 1;
    double E_0 = std::get<0>(GetParam());
    TimeStep timeData(nx);
    timeData.TotalTime = 0;
    FixedData fixedData;
    fixedData.radflxlm = 1;
    fixedData.LeftRadBoundaryCondition = "c";
    fixedData.RightRadBoundaryCondition = "c";
    MultiRadTrans radiationTransport(nx);
    radiationTransport.setPlanckianOpacity(0.04);
    radiationTransport.setRadEnergyDensity(E_0);
    std::vector<LoadInTables> opacityTable = {LoadInTables{0, fixedData.Nx}};
    opacityTable.at(0).PhotonGrid = {1.00000000e-050, 1.0e+050};
    int j = 0;

    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    ImplicitEquilibrium(&initGridData, fixedData, E_0);
    radiationTransport.storeOpacityTables(opacityTable);
    radiationTransport.setEnergyGroup(j);
    radiationTransport.planckIntegral(&initGridData, fixedData);
    radiationTransport.calcFreeFreeEmission(&initGridData, fixedData);
    radiationTransport.fullyImplicit(&initGridData, fixedData, timeData);
    radiationTransport.calcFreeFreeAbsorption(&initGridData, fixedData);
    radiationTransport.gather(&initGridData);
    timeData.Tmax = 0.5E-7;
    timeData.Dt05 = 1E-12;
    variableData.push_back(initGridData);
    std::vector<double> analytical_form = ImplicitAnalyticEquilibriumResult(E_0, timeData.Dt05, ceil(timeData.Tmax/timeData.Dt05) + 1);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {
            break;
        }
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        variableData[lastStepDataIndex] = variableData[lastStepDataIndex - 1]; // transfers specific heat, opacities, density.
        //Units being tested
        std::fill(variableData[lastStepDataIndex].RadFFEmission.begin(), variableData[lastStepDataIndex].RadFFEmission.end(), 0.0);
        std::fill(variableData[lastStepDataIndex].RadFFAbsorb.begin(), variableData[lastStepDataIndex].RadFFAbsorb.end(), 0.0);
        std::fill(variableData[lastStepDataIndex].RadEnergyDensity.begin(), variableData[lastStepDataIndex].RadEnergyDensity.end(), 0.0);
        radiationTransport.planckIntegral(&variableData[lastStepDataIndex], fixedData);
        // radiationTransport.updateRadiationEnergyDensity(&variableData[lastStepDataIndex], fixedData, timeData);
        radiationTransport.calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData);
        radiationTransport.calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.gather(&variableData.at(lastStepDataIndex));
        
        timeData.incrementStep();
        timeData.setTotalTime();
        double relative_error;
        relative_error = abs((variableData[lastStepDataIndex].RadEnergyDensity[0] - analytical_form[timeData.Step])/analytical_form[timeData.Step]);
        ASSERT_LT(relative_error, 1e-2);
    }
}

void ImplicitwriteOut(std::ofstream &variable, std::vector<double> array, int length)
{
    if (!variable)
    {
        std::cout << "file could not be open for writing ! \n";
    }
    for(int i =0; i < length; i++)
    {
        //variable << std::fixed <<  array[i] << "\n";
        variable << std::fixed << std::setprecision(35) <<array[i] << "\n";
    }
    variable.close(); 
}
std::string ImplicitfileNameCreator(std::string path, std::string var, std::string step, std::string extension)
{
    std::string updated_path = path + "/" + var + "_" + step + extension;
//    V path.append("/");
//     path.append(var);
//     path.append("_");
//     path.append(step);
//     path.append(extension);
    return(updated_path);
}
std::string ImplicitfileNameCreator(std::string path, std::string var, std::string extension)
{
    std::string updated_path = path + "/" + var + extension;
//    V path.append("/");
//     path.append(var);
//     path.append("_");
//     path.append(step);
//     path.append(extension);
    return(updated_path);
}
void ImplicitloadSuOlsenExpected(std::string path, std::vector<double> &SuOlsenU0_100,std::vector<double> &SuOlsenU0_3300,std::vector<double> &SuOlsenU1_100,
                        std::vector<double> &SuOlsenU1_3300,std::vector<double> &SuOlsenV_100,std::vector<double> &SuOlsenV_3300)
{
    std::vector<std::string> vars = {"RAD_ENERGY_0_100", "RAD_ENERGY_1_100", "ELECTRON_TEMPERATURE_100","RAD_ENERGY_0_3300", "RAD_ENERGY_1_3300",  "ELECTRON_TEMPERATURE_3300"};
    int no_var = vars.size();
    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::string file_path = ImplicitfileNameCreator(path, vars[i], ".txt");
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    SuOlsenU0_100.push_back(stod(strInput));
                }
                if(i == 1)
                {

                    SuOlsenU1_100.push_back(stod(strInput));
                }
                if(i == 2)
                {

                    SuOlsenV_100.push_back(stod(strInput));
                }
                if(i == 3)
                {

                    SuOlsenU0_3300.push_back(stod(strInput));
                }
                if(i == 4)
                {

                    SuOlsenU1_3300.push_back(stod(strInput));
                }
                if(i == 5)
                {

                    SuOlsenV_3300.push_back(stod(strInput));
                }
            }
        }
    }
}
TEST(ImplicitMultiGroupRadiationTransport, SuOlsen)
{
    int nx = 1024;
    Init init;
    TimeStep timeData(nx);
    timeData.TotalTime = 0;
    FixedData fixedData;
    fixedData.radflxlm = 0;
    fixedData.larsen_limiter = 1;
    fixedData.RightRadBoundaryCondition = "d";
    fixedData.LeftRadBoundaryCondition = "r";
    fixedData.Nx = nx;
    double kappa_1 = 2.0/101.0 * 1.0;
    double kappa_2 = 200.0/101.0 * 1.0;
    double kappa_bar = 0.5 * kappa_1+ 0.5 * kappa_2;
    
    double alpha = 4.0;
    double alpha_ratio = 4.0/alpha;
    double tao_0 = 10.0 / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    double tau = 10.0;
    double x0 = 0.5 / kappa_bar;

    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/ImplicitMGSuOlsenVari";
    std::vector<double> SuOlsenU0_100,SuOlsenU0_3300,SuOlsenU1_100,SuOlsenU1_3300,SuOlsenV_100,SuOlsenV_3300;
    ImplicitloadSuOlsenExpected(localPath, SuOlsenU0_100, SuOlsenU0_3300, SuOlsenU1_100, SuOlsenU1_3300, SuOlsenV_100, SuOlsenV_3300);

    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    InitImplicitMultiSuOlsen(&initGridData, fixedData, nx, kappa_bar, alpha, 1.0);
    // MatrixSolver matrixSolver0(nx);
    // MatrixSolver matrixSolver1(nx);
    // matrixSolver0.createPETScObjects();
    // matrixSolver1.createPETScObjects();
    Hydro.updateDxs(&initGridData);
    variableData.push_back(initGridData);
    fixedData.RadNg = 2;
    timeData.Tmax = (tau) / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    
    double timestep = (30.0) / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    timeData.Dt05 = timestep * 1e-4;
    std::vector<MultiRadTrans> radiationTransport = {MultiRadTrans(fixedData.Nx),MultiRadTrans(fixedData.Nx)};
    radiationTransport.at(0).setTotalEnergyGroup(2);
    radiationTransport.at(1).setTotalEnergyGroup(2);
    radiationTransport.at(0).setEnergyGroup(0);
    radiationTransport.at(1).setEnergyGroup(1);

    radiationTransport.at(0).setRossOpacity(kappa_1);
    radiationTransport.at(0).setPlanckianOpacity(kappa_1);
    radiationTransport.at(1).setRossOpacity(kappa_2);
    radiationTransport.at(1).setPlanckianOpacity(kappa_2);
    std::vector<double> heating;
    heating.resize(fixedData.Nx);
    double heating_temperature = 1.0e4;
    for(int i = 0; i < fixedData.Nx; i++)
    {
        if(initGridData.CellCenteredCoord[i] < x0)
        {
            heating[i] = 0.5 * kappa_bar * fixedData.SPEED_OF_LIGHT*fixedData.radiationConstant * pow(heating_temperature, 4); //kappa= 1e-2
        }
        else
        {
            heating[i] = 0.0;
        }
    }
    std::vector<double> zeroth = ImplicitMultiSuOlsenPlanckIntegral(&initGridData,fixedData);
    radiationTransport.at(0).setPlanckIntegral(zeroth);
    radiationTransport.at(1).setPlanckIntegral(zeroth);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {
            break;
        }
        if(variableData.size() > 2)
        {
            variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        //Units being tested

        variableData.at(lastStepDataIndex).Mass = variableData.at(lastStepDataIndex - 1).Mass;
        variableData.at(lastStepDataIndex).CellWallCoord = variableData.at(lastStepDataIndex - 1).CellWallCoord;
        variableData.at(lastStepDataIndex).CellCenteredCoord = variableData.at(lastStepDataIndex - 1).CellCenteredCoord;
        variableData.at(lastStepDataIndex).Density = variableData.at(lastStepDataIndex - 1).Density;
        variableData.at(lastStepDataIndex).NumberDensityE = variableData.at(lastStepDataIndex - 1).NumberDensityE;
        variableData.at(lastStepDataIndex).Dx_k_1_2 = variableData.at(lastStepDataIndex - 1).Dx_k_1_2;
        variableData.at(lastStepDataIndex).Dx_k = variableData.at(lastStepDataIndex - 1).Dx_k;
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex- 1), &variableData.at(lastStepDataIndex), timeData, false);

        std::vector<double> zeroth = ImplicitMultiSuOlsenPlanckIntegral(&variableData.at(lastStepDataIndex),fixedData);
        radiationTransport.at(0).setPlanckIntegral(zeroth);
        radiationTransport.at(1).setPlanckIntegral(zeroth);
        radiationTransport.at(0).calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).ficksDiffusionCoefficient(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).ficksDiffusionCoefficient(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData, heating);
        radiationTransport.at(1).fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData, heating);
        radiationTransport.at(0).calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).gather(&variableData.at(lastStepDataIndex));
        radiationTransport.at(1).gather(&variableData.at(lastStepDataIndex));
        
        for(int i = 0; i < nx; i++)
        {
            variableData.at(lastStepDataIndex).SpecificHeatE[i] = (alpha/variableData.at(lastStepDataIndex).Density[i])
                                                                  * fixedData.radiationConstant
                                                                  * pow(variableData.at(lastStepDataIndex).TemperatureE[i], 3);
        }

        timeData.incrementStep();
        timeData.setTotalTime();

        if((timeData.Step == 100) || (timeData.Step == 3300))
        {

            std::vector<double> relative_errors_u0;
            std::vector<double> relative_errors_u1;
            std::vector<double> relative_errors_v;
            for(int i = 0; i < fixedData.Nx; i++)
            {
                std::vector<double> U0 = radiationTransport.at(0).getLocalEnergyDensity();
                std::vector<double> U1 = radiationTransport.at(1).getLocalEnergyDensity();
                std::vector<double> V = variableData.at(lastStepDataIndex).TemperatureE;
              
                if(timeData.Step == 100)
                {
                    relative_errors_u0.push_back(abs(U0[i] - SuOlsenU0_100[i]) / SuOlsenU0_100[i]);
                    relative_errors_u1.push_back(abs(U1[i] - SuOlsenU1_100[i]) / SuOlsenU1_100[i]);
                    relative_errors_v.push_back(abs(V[i] - SuOlsenV_100[i]) / SuOlsenV_100[i]);
                }
                if(timeData.Step == 3300)
                {
                    relative_errors_u0.push_back(abs(U0[i] - SuOlsenU0_3300[i]) / SuOlsenU0_3300[i]);
                    relative_errors_u1.push_back(abs(U1[i] - SuOlsenU1_3300[i]) / SuOlsenU1_3300[i]);
                    relative_errors_v.push_back(abs(V[i] - SuOlsenV_3300[i]) / SuOlsenV_3300[i]);
                }
            }
            auto max_error1 = std::max_element(std::begin(relative_errors_u0), std::end(relative_errors_u0));
            auto max_error2 = std::max_element(std::begin(relative_errors_u1), std::end(relative_errors_u1));
            auto max_error3 = std::max_element(std::begin(relative_errors_v), std::end(relative_errors_v));
            ASSERT_LT(*max_error1, 1.0e-8);
            ASSERT_LT(*max_error2, 1.0e-8);
            ASSERT_LT(*max_error3, 1.0e-8);
        }
    }
    // matrixSolver0.cleanUpPETScObjects();
    // matrixSolver1.cleanUpPETScObjects();
}
std::vector<double> ImplicitAnalyticalMGEquilibriumTest(std::vector<double> planck, double opacity, double time, double speed_of_light)
{
    int no_photons = planck.size();
    std::vector<double> analytic_result;
    analytic_result.resize(no_photons, 0);
    for(int i = 0; i < no_photons; i++)
    {
        double rescale_planck = (planck[i] * 4 * M_PI) / speed_of_light;
        analytic_result[i] = rescale_planck * ( 1 - exp(-1* speed_of_light* opacity * time))    ;
    }
    return analytic_result;
}
TEST(ImplicitMultiGroupRadiationTransport, EquilibriumTest)
{
    FixedData fixedData;
    fixedData.Nx = 1;
    TimeStep timeData(fixedData.Nx);
    timeData.TotalTime = 0;
    fixedData.radflxlm = 1;
    fixedData.LeftRadBoundaryCondition = "c";
    fixedData.RightRadBoundaryCondition = "c";
    fixedData.RadNg = 80;

    std::vector<MultiRadTrans> radiationTransport(80, MultiRadTrans(fixedData.Nx));
    std::vector<LoadInTables> opacityTable = {LoadInTables{0, fixedData.Nx}};
    opacityTable.at(0).PhotonGrid = {1.00000000e-01, 1.16554150e-01, 1.35848699e-01, 1.58337296e-01,
       1.84548689e-01, 2.15099156e-01, 2.50706993e-01, 2.92209404e-01,
       3.40582187e-01, 3.96962672e-01, 4.62676468e-01, 5.39268624e-01,
       6.28539961e-01, 7.32589408e-01, 8.53863357e-01, 9.95213177e-01,
       1.15996226e+00, 1.35198415e+00, 1.57579363e+00, 1.83665287e+00,
       2.14069514e+00, 2.49506903e+00, 2.90810649e+00, 3.38951880e+00,
       3.95062482e+00, 4.60461718e+00, 5.36687241e+00, 6.25531252e+00,
       7.29082633e+00, 8.49776065e+00, 9.90449268e+00, 1.15440973e+01,
       1.34551244e+01, 1.56825059e+01, 1.82786114e+01, 2.13044801e+01,
       2.48312557e+01, 2.89418590e+01, 3.37329378e+01, 3.93171388e+01,
       4.58257569e+01, 5.34118215e+01, 6.22536944e+01, 7.25592643e+01,
       8.45708337e+01, 9.85708163e+01, 1.14888377e+02, 1.33907171e+02,
       1.56074365e+02, 1.81911149e+02, 2.12024994e+02, 2.47123929e+02,
       2.88033195e+02, 3.35714642e+02, 3.91289347e+02, 4.56063972e+02,
       5.31561485e+02, 6.19556971e+02, 7.22119360e+02, 8.41660082e+02,
       9.80989753e+02, 1.14338427e+03, 1.33266181e+03, 1.55327265e+03,
       1.81040373e+03, 2.11010068e+03, 2.45940991e+03, 2.86654431e+03,
       3.34107635e+03, 3.89416314e+03, 4.53880874e+03, 5.29016995e+03,
       6.16591261e+03, 7.18662703e+03, 8.37631204e+03, 9.76293929e+03,
       1.13791109e+04, 1.32628260e+04, 1.54583741e+04, 1.80173765e+04,
       2.10000000e+04};
    int j = 0;
    for(auto &i:radiationTransport)
    {
        i.storeOpacityTables(opacityTable);
        i.setRossOpacity(1);
        i.setPlanckianOpacity(1);
        i.setEnergyGroup(j);
        i.setTotalEnergyGroup(fixedData.RadNg);
        j++;
    }
    // MatrixSolver matrixSolver(fixedData.Nx);
    FluidDynamics Hydro(fixedData.Nx);
    std::vector<GridData> variableData;
    GridData initGridData(fixedData.Nx, true);
    ImplicitEquilibrium(&initGridData, fixedData, 1e-30);
    Hydro.updateDxs(&initGridData);
    initGridData.TemperatureE[0] = 1000 * fixedData.ev_to_k;
    initGridData.Density[0] = 1;
    timeData.Tmax = 0.5E-8;
    timeData.Dt05 = 1E-12;
    Hydro.updateDxs(&initGridData);
    variableData.push_back(initGridData);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {

            // variableData.back().dump(ioData, timeData.Step);
            // timeData.dump(ioData);
            break;
        }
        GridData grid(fixedData.Nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        variableData[lastStepDataIndex] = variableData[lastStepDataIndex - 1]; // transfers specific heat, opacities, density.
        //Units being tested
        for(auto &i:radiationTransport)
        {
            i.planckIntegral(&variableData.at(lastStepDataIndex),fixedData);
            i.calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
            i.fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData);
            // i.updateRadiationEnergyDensity(&variableData.at(lastStepDataIndex),fixedData, timeData);
            i.calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
            i.gather(&variableData.at(lastStepDataIndex));
        }

        timeData.incrementStep();
        timeData.setTotalTime();
        std::vector<double> all_energyDensities;
        std::vector<double> all_planck;
        if(timeData.Step % 1000 == 0)
        {

            // variableData.back().dump(ioData, timeData.Step);
            // timeData.dump(ioData);
            for(int i = 0; i < fixedData.RadNg; i++)
            {
                std::vector<double> er_local = radiationTransport.at(i).getLocalEnergyDensity();
                std::vector<double> planck_local = radiationTransport.at(i).getPlanckIntegral();
                all_energyDensities.insert(all_energyDensities.end(), er_local.begin(), er_local.end());
                all_planck.insert(all_planck.end(), planck_local.begin(), planck_local.end());
            }

            std::vector<double> expected_value = ImplicitAnalyticalMGEquilibriumTest(all_planck, 1, timeData.TotalTime, fixedData.SPEED_OF_LIGHT);
            std::vector<double> relative_errors;
            for(int i = 0; i < fixedData.RadNg; i++)
            {
                relative_errors.push_back(abs(all_energyDensities[i] - expected_value[i]) / expected_value[i]);
            }
            auto max_error = std::max_element(std::begin(relative_errors), std::end(relative_errors));
            ASSERT_LT(*max_error, 1.5e-3);
        }
    }
}

void ImplicitloadSuOlsenExpected(std::string path, std::vector<double> &SuOlsenU_100,std::vector<double> &SuOlsenU_3300, std::vector<double> &SuOlsenV_100,std::vector<double> &SuOlsenV_3300)
{
    std::vector<std::string> vars = {"RAD_ENERGY_DENSITY_100",  "ELECTRON_TEMPERATURE_100","RAD_ENERGY_DENSITY_3163",  "ELECTRON_TEMPERATURE_3163"};
    int no_var = vars.size();
    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::string file_path = ImplicitfileNameCreator(path, vars[i], ".txt");
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    SuOlsenU_100.push_back(stod(strInput));
                }
                if(i == 1)
                {

                    SuOlsenV_100.push_back(stod(strInput));
                }
                if(i == 2)
                {

                    SuOlsenU_3300.push_back(stod(strInput));
                }
                if(i == 3)
                {

                    SuOlsenV_3300.push_back(stod(strInput));
                }
            }
        }
    }
}

TEST(ImplicitGrayRadiationTransport, SuOlsen)
{
    int nx = 1024;
    Init init;
    TimeStep timeData(nx);
    timeData.TotalTime = 0;
    FixedData fixedData;
    fixedData.radflxlm = 0;
    fixedData.larsen_limiter = 1;
    fixedData.RightRadBoundaryCondition = "d";
    fixedData.LeftRadBoundaryCondition = "r";
    fixedData.Nx = nx;
    double kappa_1 = 101.0/101.0;
    double kappa_2 = 101.0/101.0;
    double kappa_bar = 0.5 * kappa_1+ 0.5 * kappa_2;
    
    double alpha = 4.0;
    double alpha_ratio = 4.0/alpha;
    double tao_0 = 10.0 / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    double tau = 10.0;
    double x0 = 0.5 / kappa_bar;

    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/ImplicitMGSuOlsenGray";
    std::vector<double> SuOlsenU_100,SuOlsenU_3163,SuOlsenV_100,SuOlsenV_3163;
    ImplicitloadSuOlsenExpected(localPath, SuOlsenU_100, SuOlsenU_3163, SuOlsenV_100, SuOlsenV_3163);

    FluidDynamics Hydro(nx);
    std::vector<GridData> variableData;
    GridData initGridData(nx, true);
    InitImplicitMultiSuOlsen(&initGridData, fixedData, nx, kappa_bar, alpha, 1.0);
    Hydro.updateDxs(&initGridData);
    fixedData.RadNg = 2;
    timeData.Tmax = (tau) / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    
    double timestep = (0.31623) / (alpha_ratio * fixedData.SPEED_OF_LIGHT * kappa_bar);
    timeData.Dt05 = timestep * 1e-2;
    std::vector<MultiRadTrans> radiationTransport = {MultiRadTrans(fixedData.Nx),MultiRadTrans(fixedData.Nx)};
    radiationTransport.at(0).setTotalEnergyGroup(2);
    radiationTransport.at(1).setTotalEnergyGroup(2);
    radiationTransport.at(0).setEnergyGroup(0);
    radiationTransport.at(1).setEnergyGroup(1);

    radiationTransport.at(0).setRossOpacity(kappa_1);
    radiationTransport.at(0).setPlanckianOpacity(kappa_1);
    radiationTransport.at(1).setRossOpacity(kappa_2);
    radiationTransport.at(1).setPlanckianOpacity(kappa_2);
    std::vector<double> heating;
    heating.resize(fixedData.Nx);
    double heating_temperature = 1.0e4;
    for(int i = 0; i < fixedData.Nx; i++)
    {
        if(initGridData.CellCenteredCoord[i] < x0)
        {
            heating[i] = 0.5 * kappa_bar * fixedData.SPEED_OF_LIGHT*fixedData.radiationConstant * pow(heating_temperature, 4); //kappa= 1e-2
        }
        else
        {
            heating[i] = 0.0;
        }
    }
    std::vector<double> zeroth = ImplicitMultiSuOlsenPlanckIntegral(&initGridData,fixedData);
    radiationTransport.at(0).setPlanckIntegral(zeroth);
    radiationTransport.at(1).setPlanckIntegral(zeroth);
    variableData.push_back(initGridData);
    while(true)
    {
        if(timeData.TotalTime >= timeData.Tmax)
        {

            break;
        }
        if(variableData.size() > 2)
        {
            variableData.erase(variableData.begin());   
        } 
        GridData grid(nx);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        //Units being tested

        variableData.at(lastStepDataIndex).Mass = variableData.at(lastStepDataIndex - 1).Mass;
        variableData.at(lastStepDataIndex).CellWallCoord = variableData.at(lastStepDataIndex - 1).CellWallCoord;
        variableData.at(lastStepDataIndex).CellCenteredCoord = variableData.at(lastStepDataIndex - 1).CellCenteredCoord;
        variableData.at(lastStepDataIndex).Density = variableData.at(lastStepDataIndex - 1).Density;
        variableData.at(lastStepDataIndex).NumberDensityE = variableData.at(lastStepDataIndex - 1).NumberDensityE;
        variableData.at(lastStepDataIndex).Dx_k_1_2 = variableData.at(lastStepDataIndex - 1).Dx_k_1_2;
        variableData.at(lastStepDataIndex).Dx_k = variableData.at(lastStepDataIndex - 1).Dx_k;
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex- 1), &variableData.at(lastStepDataIndex), timeData, false);

        std::vector<double> zeroth = ImplicitMultiSuOlsenPlanckIntegral(&variableData.at(lastStepDataIndex),fixedData);
        radiationTransport.at(0).setPlanckIntegral(zeroth);
        radiationTransport.at(1).setPlanckIntegral(zeroth);
        radiationTransport.at(0).calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).calcFreeFreeEmission(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).ficksDiffusionCoefficient(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).ficksDiffusionCoefficient(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData, heating);
        radiationTransport.at(1).fullyImplicit(&variableData.at(lastStepDataIndex), fixedData, timeData, heating);
        radiationTransport.at(0).calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(1).calcFreeFreeAbsorption(&variableData.at(lastStepDataIndex), fixedData);
        radiationTransport.at(0).gather(&variableData.at(lastStepDataIndex));
        radiationTransport.at(1).gather(&variableData.at(lastStepDataIndex));
        
        for(int i = 0; i < nx; i++)
        {
            variableData.at(lastStepDataIndex).SpecificHeatE[i] = (alpha/variableData.at(lastStepDataIndex).Density[i])
                                                                  * fixedData.radiationConstant
                                                                  * pow(variableData.at(lastStepDataIndex).TemperatureE[i], 3);
        }

        timeData.incrementStep();
        timeData.setTotalTime();

        if((timeData.Step == 100) || (timeData.Step == 3163))
        {

            std::vector<double> relative_errors_u0;
            std::vector<double> relative_errors_u1;
            std::vector<double> relative_errors_v;

            for(int i = 0; i < fixedData.Nx; i++)
            {
                std::vector<double> U = variableData.at(lastStepDataIndex).RadEnergyDensity; 
                std::vector<double> V = variableData.at(lastStepDataIndex).TemperatureE;
                SuOlsenU_100[fixedData.Nx - 1] = 1e-15; // accounts for the changes to floor val  
                SuOlsenU_3163[fixedData.Nx - 1] = 1e-15;  
                if(timeData.Step == 100)
                {
                    relative_errors_u1.push_back(abs(U[i] - SuOlsenU_100[i]) / SuOlsenU_100[i]);
                    relative_errors_v.push_back(abs(V[i] - SuOlsenV_100[i]) / SuOlsenV_100[i]);
                }
                if(timeData.Step == 3163)
                {
                    relative_errors_u1.push_back(abs(U[i] - SuOlsenU_3163[i]) / SuOlsenU_3163[i]);
                    relative_errors_v.push_back(abs(V[i] - SuOlsenV_3163[i]) / SuOlsenV_3163[i]);
                }
            }
            auto max_error2 = std::max_element(std::begin(relative_errors_u1), std::end(relative_errors_u1));
            auto max_error3 = std::max_element(std::begin(relative_errors_v), std::end(relative_errors_v));
            ASSERT_LT(*max_error2, 1.0e-8);
            ASSERT_LT(*max_error3, 1.0e-8);
        }
    }
}
INSTANTIATE_TEST_SUITE_P(
    ImplicitGrayRadiationTransportEquili,
    ImplicitGrayRadiationTransportEquil,
    ::testing::Values(
        // std::make_tuple(1E7),
        std::make_tuple(1E11)));
