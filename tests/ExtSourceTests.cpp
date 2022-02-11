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
#include "gtest/gtest.h"
#include "HyKiCT/FixedData.h"
#include "HyKiCT/GridData.h"
#include "HyKiCT/ExtSource.h"
#include "HyKiCT/FluidDynamics.h"
#include "HyKiCT/TimeStep.h"
#include "HyKiCT/IdealGasEoS.h"
#include "HyKiCT/PlasmaParameters.h"
#include <iostream> 
#include <numeric>
#include <algorithm>

class VariableLaserPowerChecks:public::testing::TestWithParam<std::tuple<std::string>>
{
};

void initEnergyConsWithCrit(GridData *gridData, FixedData const &fixedData)
{
    double dx = 0.01;
    int nx = 100;
    for(int i =0; i < nx + 1; i++)
    {
        gridData->CellWallCoord[i] = i*dx;
    }
    for(int i = 0; i < nx; i++)
    {
        gridData->CellCenteredCoord[i] = (gridData->CellWallCoord[i + 1] + gridData->CellWallCoord[i]) / 2;
        gridData->NumberDensityI[i] = gridData->NumberDensityE[i] / gridData->Zbar[i];
        gridData->Density[i] = gridData->NumberDensityI[i] * fixedData.Ar[i] * fixedData.PROTON_MASS;
        gridData->Mass[i] = gridData->Density[i] * (gridData->CellWallCoord[i - 1] - gridData->CellWallCoord[i]);
    }
}
void initLaserPropagation(GridData *gridData, FixedData const &fixedData)
{
    int nx = 100;
    double dx = 0.30e-2/(nx+1);
    for(int i =0; i < nx + 1; i++)
    {
        gridData->CellWallCoord[i] = i*dx;
    }
    for(int i = 0; i < nx; i++)
    {
        gridData->CellCenteredCoord[i] = (gridData->CellWallCoord[i + 1] + gridData->CellWallCoord[i]) / 2;
    
        // if(gridData->CellCenteredCoord[i] > 0.3e-2)
        // {
        //     gridData->Density[i] *=0.00001;
        // }
        gridData->Mass[i] = gridData->Density[i] * (gridData->CellWallCoord[i + 1] - gridData->CellWallCoord[i]);
    }
}
TEST(LaserPowerConsTests, EnergyConservationWithCritical)
{
    int nx = 100;
    int critcalSurface = 24; // std::get<0>(GetParam());
    FixedData fixedData;
    GridData TestGridData(nx);
    ExtSource laser;
    laser.setStartIndex(nx); 
    laser.mNx = nx;
    PlasmaParameters plasmaParams(nx);
    laser.setLaserPower(1e10);
    std::fill(TestGridData.Zbar.begin(), TestGridData.Zbar.end(), 30);
    std::fill(TestGridData.NumberDensityE.begin(), TestGridData.NumberDensityE.end(), 1e24);
    std::fill(TestGridData.TemperatureE.begin(), TestGridData.TemperatureE.end(), 100);
    std::fill(TestGridData.CoulombLogEI.begin(), TestGridData.CoulombLogEI.end(), 2.5);
    fixedData.Nx = 100;
    fixedData.Ar.resize(nx, 59.9);
    laser.setLaserWavelength(100e-9);
    initEnergyConsWithCrit(&TestGridData, fixedData);
    TestGridData.NumberDensityE[critcalSurface] = 1e31;
    plasmaParams.calculateCoulombLogLaser(&TestGridData, fixedData, laser.LaserWavelength);
    plasmaParams.calculateCollisionFrequencyEIOverC(&TestGridData, fixedData, true);
    laser.calculateBeta(&TestGridData);
    // laser.inverseBrem(&TestGridData, fixedData);
    laser.inverseBrem(&TestGridData);
    double PowerSum = 0;  // 0 - value of 3rd param
    for (auto x : TestGridData.PowerAbsorbed) PowerSum += x;
    double PowerLoss = abs(PowerSum - 1e10) / 1e10;
    ASSERT_LT(PowerLoss, 1e-5);
}
TEST(LaserPowerConsTests, VarLaserPowerEnergyConservationWithCritical)
{
    std::cout << "\n New Test" << "\n";
    int nx = 100;
    int nt = 100;
    double dt = 1e-15;
    int critcalSurface = 24; //std::get<0>(GetParam());
    FixedData fixedData;
    GridData TestGridData(nx);
    PlasmaParameters plasmaParams(nx);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LOAD_LASER_TESTS/exp_test_laser.txt";
    ExtSource laser;
    laser.setStartIndex(nx); 
    laser.Init(localPath);
    laser.mNx = nx;
    std::fill(TestGridData.Zbar.begin(), TestGridData.Zbar.end(), 30);
    std::fill(TestGridData.NumberDensityE.begin(), TestGridData.NumberDensityE.end(), 1e24);
    std::fill(TestGridData.TemperatureE.begin(), TestGridData.TemperatureE.end(), 100);
    std::fill(TestGridData.CoulombLogEI.begin(), TestGridData.CoulombLogEI.end(), 2.5);
    fixedData.Nx = 100;
    fixedData.Ar.resize(nx, 59.9);
    laser.setLaserWavelength(100e-9);
    initEnergyConsWithCrit(&TestGridData, fixedData);
    TestGridData.NumberDensityE[critcalSurface] = 1e31;
    plasmaParams.calculateCoulombLogLaser(&TestGridData, fixedData, laser.LaserWavelength);
    plasmaParams.calculateCollisionFrequencyEIOverC(&TestGridData, fixedData, true);
    laser.calculateBeta(&TestGridData);
    for(int i = 0; i < nt; i++ )
    {
        laser.updateLaserPower(i * dt);
        double object_laser= laser.getLaserPower();
        // laser.inverseBrem(&TestGridData, fixedData);
        if(object_laser > 0)
        {
            laser.inverseBrem(&TestGridData);
            double PowerSum = 0;  // 0 - value of 3rd param
            for (auto x : TestGridData.PowerAbsorbed) PowerSum += x;
            double PowerLoss = abs(PowerSum - object_laser) / object_laser;
            std::fill(TestGridData.PowerAbsorbed.begin(), TestGridData.PowerAbsorbed.end(), 0);
            ASSERT_LT(PowerLoss, 1e-4);
        }
    }
}


TEST(LoadCheck, VariableLoadLaser)
{
    int nx = 100;
    int nt = 900;
    double dt = 1e-16;
    std::string profile_path = "exp_test_laser.txt";
    FixedData fixedData;
    GridData TestGridData(nx);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LOAD_LASER_TESTS/" + profile_path;
    ExtSource laser;
    laser.Init(localPath);
    laser.mNx = nx;
    laser.setStartIndex(nx); 
    std::vector<double> true_laser;
    if(profile_path == "exp_test_laser.txt")
    {
        true_laser = {
                    0, 9.048374180359594345E9, 8.187307530779818535E9, 7.408182206817178726E9, 6.703200460356392860E9,
                    6.065306597126333237E9, 5.488116360940263748E9, 4.965853037914094925E9, 4.493289641172215462E9, 4.065696597405990601E9, 3.678794411714422226E9, 3.328710836980795383E9,
                    3.011942119122020245E9, 2.725317930340126038E9, 2.465969639416064262E9, 2.231301601484297752E9, 2.018965179946553946E9, 1.826835240527346611E9, 1.652988882215865135E9,
                    1.495686192226350307E9, 1.353352832366126299E9, 1.224564282529819012E9, 1.108031583623338699E9, 1.002588437228037119E9, 9.071795328941247463E8, 8.208499862389880419E8,
                    7.427357821433387995E8, 6.720551273974975348E8, 6.081006262521795034E8, 5.502322005640721321E8, 4.978706836786392331E8, 4.504920239355779886E8, 4.076220397836621404E8,
                    3.688316740123999119E8, 3.337326996032608151E8, 3.019738342231848836E8, 2.732372244729254246E8, 2.472352647033938766E8, 2.237077185616559088E8, 2.024191144580438137E8,
                    1.831563888873416483E8, 1.657267540176123679E8, 1.499557682047770321E8, 1.356855901220092177E8, 1.227733990306843668E8, 1.110899653824230582E8, 1.005183574463357627E8,
                    9.095277101695816219E7, 8.229747049020023644E7, 7.446583070924338698E7, 6.737946999085466564E7, 6.096746565515632927E7, 5.516564420760771632E7, 4.991593906910216808E7,
                    4.516580942612665892E7, 4.086771438464063406E7, 3.697863716482929140E7, 3.345965457471268997E7, 3.027554745375812799E7, 2.739444818768368289E7, 2.478752176666356251E7,
                    2.242867719485801086E7, 2.029430636295733973E7, 1.836304777028905600E7, 1.661557273173933849E7, 1.503439192977572419E7, 1.360368037547892705E7, 1.230911902673481032E7, 
                    1.113775147844803333E7, 1.007785429048509523E7, 9.118819655545154586E6, 8.251049232659045607E6, 7.465858083766785450E6, 6.755387751938437112E6, 6.112527611295723356E6, 
                    5.530843701478330418E6, 5.004514334406103939E6, 4.528271828867969103E6, 4.097349789797864389E6, 3.707435404590882361E6, 3.354626279025112744E6, 3.035391380788662471E6, 
                    2.746535699721420649E6, 2.485168271079518367E6, 2.248673241788481828E6, 2.034683690106441732E6, 1.841057936675788602E6, 1.665858109876332339E6, 1.507330750954765128E6,
                    1.363889264820114011E6, 1.234098040866795694E6, 1.116658084901147755E6, 1.010394018370932434E6, 9.142423147817327408E5, 8.272406555663222680E5, 7.485182988770046504E5,
                    6.772873649085377110E5, 6.128349505322201876E5, 5.545159943217694527E5, 5.017468205617528292E5};
    }
    else
    {
        true_laser = {
            0, 6.279051952931337357e+08, 1.253332335643042564e+09, 1.873813145857246399e+09, 2.486898871648548126e+09, 3.090169943749474525e+09, 
            3.681245526846779823e+09, 4.257792915650725842e+09, 4.817536741017152786e+09, 5.358267949789966583e+09, 5.877852522924732208e+09, 6.374239897486896515e+09, 
            6.845471059286887169e+09, 7.289686274214115143e+09, 7.705132427757891655e+09, 8.090169943749474525e+09, 8.443279255020151138e+09, 8.763066800438634872e+09, 
            9.048270524660196304e+09, 9.297764858882513046e+09, 9.510565162951536179e+09, 9.685831611286310196e+09, 9.822872507286886215e+09, 9.921147013144779205e+09, 
            9.980267284282714844e+09, 1.000000000000000000e+10, 9.980267284282714844e+09, 9.921147013144777298e+09, 9.822872507286888123e+09, 9.685831611286310196e+09, 
            9.510565162951536179e+09, 9.297764858882513046e+09, 9.048270524660194397e+09, 8.763066800438634872e+09, 8.443279255020152092e+09, 8.090169943749474525e+09, 
            7.705132427757892609e+09, 7.289686274214115143e+09, 6.845471059286888123e+09, 6.374239897486898422e+09, 5.877852522924728394e+09, 5.358267949789966583e+09, 
            4.817536741017151833e+09, 4.257792915650724888e+09, 3.681245526846781254e+09, 3.090169943749475002e+09, 2.486898871648548126e+09, 1.873813145857245445e+09, 
            1.253332335643040895e+09, 6.279051952931313515e+08, 1.224646799147353167e-06, 6.279051952931333780e+08, 1.253332335643043041e+09, 1.873813145857247591e+09, 
            2.486898871648550034e+09, 3.090169943749477386e+09, 3.681245526846775532e+09, 4.257792915650727272e+09, 4.817536741017153740e+09, 5.358267949789964676e+09, 
            5.877852522924733162e+09, 6.374239897486896515e+09, 6.845471059286887169e+09, 7.289686274214118958e+09, 7.705132427757893562e+09, 8.090169943749473572e+09, 
            8.443279255020153046e+09, 8.763066800438636780e+09, 9.048270524660194397e+09, 9.297764858882514954e+09, 9.510565162951536179e+09, 9.685831611286310196e+09, 
            9.822872507286888123e+09, 9.921147013144777298e+09, 9.980267284282714844e+09, 1.000000000000000000e+10, 9.980267284282714844e+09, 9.921147013144779205e+09, 
            9.822872507286888123e+09, 9.685831611286310196e+09, 9.510565162951532364e+09, 9.297764858882513046e+09, 9.048270524660196304e+09, 8.763066800438632965e+09, 
            8.443279255020151138e+09, 8.090169943749475479e+09, 7.705132427757890701e+09, 7.289686274214108467e+09, 6.845471059286889076e+09, 6.374239897486896515e+09, 
            5.877852522924732208e+09, 5.358267949789962769e+09, 4.817536741017152786e+09, 4.257792915650722027e+09, 3.681245526846777916e+09, 3.090169943749475956e+09, 
            2.486898871648544788e+09, 1.873813145857246876e+09, 1.253332335643037796e+09, 6.279051952931326628e+08};
    }
    double object_laser;
    laser.updateLaserPower(1.3e-15);
    object_laser= laser.getLaserPower();
    ASSERT_LT(object_laser, true_laser[1]);
    ASSERT_LT(true_laser[0], object_laser );


    laser.updateLaserPower(5e-15 + 3e-16);
    object_laser= laser.getLaserPower();
    ASSERT_LT(object_laser, true_laser[5]);
    ASSERT_LT(true_laser[6], object_laser );

    laser.updateLaserPower(7e-14 + 5e-16);
    object_laser= laser.getLaserPower();
    ASSERT_LT(object_laser, true_laser[70]);
    ASSERT_LT(true_laser[71], object_laser );

    laser.updateLaserPower(9e-14 + 8e-15);
    object_laser= laser.getLaserPower();
    ASSERT_LT(object_laser, true_laser[97]);
    ASSERT_LT(true_laser[98], object_laser );
}
TEST_P(VariableLaserPowerChecks, EnergyConservationWithVariablePower)
{
    int nx = 100;
    int nt = 100;
    double dt = 1e-15;
    std::string profile_path = std::get<0>(GetParam());
    FixedData fixedData;
    GridData TestGridData(nx);
    PlasmaParameters plasmaParams(nx);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LOAD_LASER_TESTS/" + profile_path;
    ExtSource laser;
    laser.Init(localPath);
    laser.mNx = nx;
    laser.setStartIndex(nx); 
    std::fill(TestGridData.Zbar.begin(), TestGridData.Zbar.end(), 30);
    std::fill(TestGridData.NumberDensityE.begin(), TestGridData.NumberDensityE.end(), 1e24);
    std::fill(TestGridData.TemperatureE.begin(), TestGridData.TemperatureE.end(), 100);
    std::fill(TestGridData.CoulombLogEI.begin(), TestGridData.CoulombLogEI.end(), 2.5);
    fixedData.Nx = 100;
    fixedData.Ar.resize(nx, 59.9);
    laser.setLaserWavelength(100e-9);
    initEnergyConsWithCrit(&TestGridData, fixedData);
    TestGridData.NumberDensityE[0] = 1e31;
    plasmaParams.calculateCoulombLogLaser(&TestGridData, fixedData, laser.LaserWavelength);
    plasmaParams.calculateCollisionFrequencyEIOverC(&TestGridData, fixedData, true);
    laser.calculateBeta(&TestGridData);
    for(int i = 0; i < nt; i++ )
    {
        laser.updateLaserPower(i * dt);
        double laserPower = laser.getLaserPower();
        // laser.inverseBrem(&TestGridData, fixedData);
        laser.inverseBrem(&TestGridData);
        double PowerSum = 0;  // 0 - value of 3rd param
        for (auto x : TestGridData.PowerAbsorbed) PowerSum += x;
        if(laserPower < 1e-5)
        {
            std::cout << PowerSum << "\n";
            std::cout << laserPower << "\n";
            continue;
        }
        double PowerLoss = abs(PowerSum - laserPower) / laserPower;
        std::fill(TestGridData.PowerAbsorbed.begin(), TestGridData.PowerAbsorbed.end(), 0);
        ASSERT_LT(PowerLoss, 1e-4);
    }
}

TEST(LaserTimeDevelop, TemporalAccuracy)
{
    int nx = 1;
    int Z = 1; //std::get<0>(GetParam());
    FixedData fixedData;
    GridData TestGridData(nx);
    ExtSource laser;
    laser.mNx = nx;
    laser.setStartIndex(nx); 
    FluidDynamics Hydro(nx);
    IdealGasEoS EoS(0, nx, fixedData.Gamma, fixedData.BOLTZMANN_CONSTANT);
    TimeStep timeData(nx);
    std::vector<GridData> variableData;
    PlasmaParameters plasmaParams(nx);
    laser.setLaserPower(1e15);
    TestGridData.CellWallCoord = {0, 1e-19};
    TestGridData.CellCenteredCoord[0] = (TestGridData.CellWallCoord[1] + TestGridData.CellWallCoord[0]) / 2;
    std::fill(TestGridData.Zbar.begin(), TestGridData.Zbar.end(), Z);
    std::fill(TestGridData.NumberDensityE.begin(), TestGridData.NumberDensityE.end(), 9.98E26);
    std::fill(TestGridData.TemperatureE.begin(), TestGridData.TemperatureE.end(), 100);
    std::fill(TestGridData.CoulombLogEI.begin(), TestGridData.CoulombLogEI.end(), 10);
    std::fill(TestGridData.Mass.begin(), TestGridData.Mass.end(), 1.67E-19);
    std::fill(TestGridData.DpDtE.begin(), TestGridData.DpDtE.end(), 0);
    std::fill(TestGridData.Density.begin(), TestGridData.Density.end(), 1.67);
    // std::fill(TestGridData.Density.begin(), TestGridData.Density.end(), 1.67);

    EoS.updateSpecificHeatCapacity(&TestGridData);
    
    fixedData.Nx = nx;
    fixedData.Ar.resize(nx, 2*Z);
    laser.setLaserWavelength(100e-9);
    double criticalDensity = 1114326918632954.5 / std::pow(100e-9, 2);
    double beta = 9.98E26 / criticalDensity;
    double kk = std::pow(beta, 2) * std::pow((1 - beta), -0.5);
    laser.setBeta(std::vector<double>{kk, kk});
    laser.setNoFullDump(true);
    plasmaParams.setCoulombLog(&TestGridData, 10);
    // plasmaParams.calculateCoulombLogLaser(&TestGridData, fixedData, laser.LaserWavelength);
    plasmaParams.calculateCollisionFrequencyEIOverC(&TestGridData, fixedData, false);
    std::fill(TestGridData.InterCoulombLog.begin(), TestGridData.InterCoulombLog.end(), 10);
    // laser.inverseBrem(&TestGridData, fixedData);
    laser.inverseBrem(&TestGridData);
    variableData.push_back(TestGridData);
    
    timeData.Dt05 = 1e-16;
    timeData.TotalTime = 0.0;
    timeData.MaxSteps = 3000;
    for(int t = 0; t < timeData.MaxSteps; t++ )
    {
        if(variableData.size() > 4)
        {
                variableData.erase(variableData.begin());   
        } 
        GridData grid(nx, false);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1;  // index values is -1 of size
        
        variableData.at(lastStepDataIndex).CellWallCoord =  variableData.at(lastStepDataIndex - 1).CellWallCoord; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).CellCenteredCoord =  variableData.at(lastStepDataIndex - 1).CellCenteredCoord; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).Mass =  variableData.at(lastStepDataIndex - 1).Mass; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).NumberDensityE =  variableData.at(lastStepDataIndex - 1).NumberDensityE; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).Zbar =  variableData.at(lastStepDataIndex - 1).Zbar; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).CoulombLogEI =  variableData.at(lastStepDataIndex - 1).CoulombLogEI; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).DpDtE =  variableData.at(lastStepDataIndex - 1).DpDtE; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).Density =  variableData.at(lastStepDataIndex - 1).Density; //Mass stays fixed between steps 
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);
        EoS.updateSpecificHeatCapacity(&variableData.at(lastStepDataIndex));
        double beta = variableData.at(lastStepDataIndex).NumberDensityE[0] / criticalDensity;
        double alpha = 1.2130370535544696E-14 * variableData.at(lastStepDataIndex).NumberDensityE[0] * variableData.at(lastStepDataIndex).Zbar[0] * variableData.at(lastStepDataIndex).CoulombLogEI[0] 
                         * std::pow(beta, 2) * std::pow((1 - beta), -0.5);
        timeData.setTotalTime();
        double analytic_sol = pow((((2.5 * 1e15 * timeData.TotalTime * alpha * 1e-19) / (variableData.at(lastStepDataIndex).SpecificHeatE[0] * variableData.at(lastStepDataIndex).Mass[0])) + 
                                pow(TestGridData.TemperatureE[0], 2.5)), 0.4);
        plasmaParams.setCoulombLog(&variableData.at(lastStepDataIndex), 10);
        // plasmaParams.calculateCoulombLogLaser(&variableData.at(lastStepDataIndex), fixedData, laser.LaserWavelength);
        plasmaParams.calculateCollisionFrequencyEIOverC(&variableData.at(lastStepDataIndex), fixedData, false);
        std::fill(variableData.at(lastStepDataIndex).InterCoulombLog.begin(), variableData.at(lastStepDataIndex).InterCoulombLog.end(), 10);
        // laser.inverseBrem(&variableData.at(lastStepDataIndex), fixedData);
        laser.inverseBrem(&variableData.at(lastStepDataIndex));
        double diff =  abs(analytic_sol - variableData.at(lastStepDataIndex).TemperatureE[0]);
        double rel_diff = diff/analytic_sol;
        timeData.incrementStep();
        ASSERT_LT(rel_diff , 0.1);
    }
}
TEST(LaserPropagation, DISABLED_TemporalAccuracy)
{
    IO ioData;  

    int nx = 100;
    int Z = 6; 
    FixedData fixedData;
    fixedData.Nx = nx;
    GridData initGridData(nx);
    ExtSource laser;
    laser.mNx = nx;
    laser.setNoFullDump(true);
    FluidDynamics Hydro(nx);
    IdealGasEoS EoS(0, nx, 5.0/3.0, fixedData.BOLTZMANN_CONSTANT);
    TimeStep timeData(nx);
    timeData.Dt05 = 1e-14;
    timeData.Dt1 = 1e-14;
    int no_time_steps = 1000000;///1.0e-9 / timeData.Dt05;
    std::vector<GridData> variableData;
    PlasmaParameters plasmaParams(nx);
    laser.setLaserWavelength(0.53e-6);
    double n_critical = 1114326918632954.5 / std::pow(0.53e-6, 2);
    // laser.setNoFullDump(true);
    std::vector<double> intensity{0.00000000e+00, 1.30612245e+16, 2.61224490e+16, 3.91836735e+16,
       5.22448980e+16, 6.53061224e+16, 7.83673469e+16, 9.14285714e+16,
       1.04489796e+17, 1.17551020e+17, 1.30612245e+17, 1.43673469e+17,
       1.56734694e+17, 1.69795918e+17, 1.82857143e+17, 1.95918367e+17,
       2.08979592e+17, 2.22040816e+17, 2.35102041e+17, 2.48163265e+17,
       2.61224490e+17, 2.74285714e+17, 2.87346939e+17, 3.00408163e+17,
       3.13469388e+17, 3.26530612e+17, 3.39591837e+17, 3.52653061e+17,
       3.65714286e+17, 3.78775510e+17, 3.91836735e+17, 4.04897959e+17,
       4.17959184e+17, 4.31020408e+17, 4.44081633e+17, 4.57142857e+17,
       4.70204082e+17, 4.83265306e+17, 4.96326531e+17, 5.09387755e+17,
       5.22448980e+17, 5.35510204e+17, 5.48571429e+17, 5.61632653e+17,
       5.74693878e+17, 5.87755102e+17, 6.00816327e+17, 6.13877551e+17,
       6.26938776e+17, 6.40000000e+17};
    std::vector<double> time(1000, 0);
    int j = 0;
    int step = no_time_steps/1000;
    for(int i = 0; i < no_time_steps; i++)
    {
        if(i % step == 0)
        {
            time[j] = i *timeData.Dt05;
            if(time[j] >= 0.05E-9)
            {
                intensity.push_back(0.64E18);
            }
            j++;
        }
    }
    laser.setLaserPowerOverTime(intensity, time);
    fixedData.Ar.resize(nx, 12);
    std::fill(initGridData.Zbar.begin(), initGridData.Zbar.end(), Z);
    std::fill(initGridData.Density.begin(), initGridData.Density.end(), 1);
    std::fill(initGridData.TemperatureE.begin(), initGridData.TemperatureE.end(), 1*fixedData.ev_to_k);
    initLaserPropagation(&initGridData, fixedData);
    Hydro.updateDxs(&initGridData);
    Hydro.updateNumberDensityI(&initGridData, fixedData);
    Hydro.updateNumberDensityE(&initGridData);
    EoS.updatePressureTerms(&initGridData);
    EoS.updateSpecificHeatCapacity(&initGridData);
    std::fill(initGridData.CoulombLogEI.begin(), initGridData.CoulombLogEI.end(), 8);
    std::vector<double> beta(100, 0.0);
    for(int i = 0; i < nx; i++)
    {   
        beta[i] = initGridData.NumberDensityE[i] / n_critical;
    }
    laser.setBeta(beta);
    laser.updateLaserPower(timeData.TotalTime);
    plasmaParams.calculatePlasmaFrequency(&initGridData, fixedData);
    plasmaParams.calculateMinImpactParameter(&initGridData, fixedData);
    plasmaParams.calculateThermalVelocity(&initGridData, fixedData);
    // plasmaParams.calculateCoulombLogLaser(&initGridData, fixedData, laser.LaserWavelength);
    plasmaParams.calculateCollisionFrequencyEIOverC(&initGridData, fixedData, true);
    laser.calculateBeta(&initGridData);
    laser.inverseBrem(&initGridData);
    variableData.push_back(initGridData);
    
    timeData.TotalTime = 0;
    timeData.MaxSteps = no_time_steps;
    for(int t = 0; t < timeData.MaxSteps; t++ )
    {
        if(variableData.size() > 4)
        {
                variableData.erase(variableData.begin());   
        } 
        GridData grid(nx, false);
        variableData.push_back(grid);

        int lastStepDataIndex = variableData.size() - 1; 
        
        variableData.at(lastStepDataIndex).Mass =  variableData.at(lastStepDataIndex - 1).Mass; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).Zbar =  variableData.at(lastStepDataIndex - 1).Zbar; //Mass stays fixed between steps 
        variableData.at(lastStepDataIndex).CoulombLogEI =  variableData.at(lastStepDataIndex - 1).CoulombLogEI; //Mass stays fixed between steps 

        // Hydro.updateViscosity(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), constantData);
        Hydro.updateMomentum(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), fixedData, timeData);
        Hydro.updateCoords(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData);
        Hydro.updateCellCentreCoords(&variableData.at(lastStepDataIndex));
        Hydro.updateDxs(&variableData.at(lastStepDataIndex));
        Hydro.updateDensity(&variableData.at(lastStepDataIndex));
        Hydro.updateNumberDensityI(&variableData.at(lastStepDataIndex), fixedData);
        Hydro.updateNumberDensityE(&variableData.at(lastStepDataIndex));
        Hydro.updateTemperatureE(&variableData.at(lastStepDataIndex - 1), &variableData.at(lastStepDataIndex), timeData, false);
        
        EoS.updatePressureTerms(&variableData.at(lastStepDataIndex));
        EoS.updateSpecificHeatCapacity(&variableData.at(lastStepDataIndex));
        laser.updateLaserPower(timeData.TotalTime);
        plasmaParams.calculatePlasmaFrequency(&variableData.at(lastStepDataIndex), fixedData);
        plasmaParams.calculateMinImpactParameter(&variableData.at(lastStepDataIndex), fixedData);
        plasmaParams.calculateThermalVelocity(&variableData.at(lastStepDataIndex), fixedData);
        // plasmaParams.calculateCoulombLogLaser(&variableData.at(lastStepDataIndex), fixedData, laser.LaserWavelength);
        plasmaParams.calculateCollisionFrequencyEIOverC(&variableData.at(lastStepDataIndex), fixedData, true);
        // laser.calculateBeta(&variableData.at(lastStepDataIndex));
        for(int i = 0; i < nx; i++)
        {   
            beta[i] = variableData.at(lastStepDataIndex).NumberDensityE[i] / n_critical;
        }
        laser.setBeta(beta);
        laser.inverseBrem(&variableData.at(lastStepDataIndex));
        timeData.incrementStep();
        timeData.setTotalTime();
    }
}
INSTANTIATE_TEST_SUITE_P(
    LaserPowerChecks,
    VariableLaserPowerChecks,
    ::testing::Values(
        std::make_tuple("exp_test_laser.txt"),
        std::make_tuple("sin_test_laser.txt")));
