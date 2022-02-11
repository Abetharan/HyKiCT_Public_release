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
#include "HyKiCT/FixedData.h"

void FixedData::initVectors(std::string varFilePath, bool coupleMulti)
{
    std::cout << varFilePath << std::endl;
    int no_var;
    std::vector<std::string> vars;
    if ((coupleMulti && (FrontHeatStartIndex > 0)) || (PreHeatStartIndex > 0))
    {
        vars = {"Z", "Ar", "pre_heat_fit_param", "front_heat_fit_param"};
        no_var = 4;
    }
    else
    {
        vars = {"Z", "Ar"};
        no_var = 2;
    }
    

    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::cout << varFilePath << '\n';
        std::cout << vars.at(i) << '\n';

        std::string file_path = varFilePath + '/'  + vars.at(i) +  ".txt";
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    Z.push_back(stod(strInput));
                }
                if (i == 1)
                {
                    Ar.push_back(stod(strInput));
                }
                if (i == 2)
                {
                    PreB.push_back(stod(strInput));
                }
                if (i == 3)
                {
                    FrontB.push_back(stod(strInput));
                }
            }
        }
    }
}
void FixedData::loadLeakAreas(std::string path)
{
    std::cout << "Loading Leak Profile" << std::endl;
    std::ifstream initFile(path);
    std::string line;

    if (initFile.is_open())
    {
        double time;
        double area1,area2;
        
        while (getline(initFile, line))
        {
            std::istringstream ss(line);
            ss >> time >> area1 >> area2;
            leakAreaMultiplier1.push_back(area1);
            leakAreaMultiplier2.push_back(area2);
            leakTimes.push_back(time);
        }
    }
}
void FixedData::init(std::string path)
{

    YAML::Node config = YAML::LoadFile(path);

    Nx = config["FixedParameters"]["nx"].as<int>();
    Ng = config["FixedParameters"]["ng"].as<int>();
    Np = config["FixedParameters"]["np"].as<int>();
    RadNg = config["FixedParameters"]["RadNg"].as<int>();
    MaxE = config["FixedParameters"]["MaxE"].as<double>();
    PreHeatStartIndex = config["FixedParameters"]["Preheat_StartIndex"].as<int>();
    PreHeatLastIndex = config["FixedParameters"]["Preheat_LastIndex"].as<int>();
    FrontHeatStartIndex = config["FixedParameters"]["Frontheat_StartIndex"].as<int>();
    FrontHeatLastIndex = config["FixedParameters"]["Frontheat_LastIndex"].as<int>();
    OutputFrequency = config["FixedParameters"]["Output_Frequency"].as<int>();
    ConstantOpacity = config["FixedParameters"]["FixedOpacityValue"].as<double>();
    heatflxlmE = config["Limiters"]["heatflxlmE"].as<float>();
    heatflxlmI = config["Limiters"]["heatflxlmI"].as<float>();
    radflxlm = config["Limiters"]["radflxlm"].as<float>();
    prrm1 = config["Limiters"]["prrm1"].as<float>();
    prrm2 = config["Limiters"]["prrm2"].as<float>();
    radabm = config["Limiters"]["radabm"].as<float>();
    radbrm = config["Limiters"]["radbrm"].as<float>();
    radcmm = config["Limiters"]["radcmm"].as<float>();
    raddfm = config["Limiters"]["raddfm"].as<float>();
    radesm = config["Limiters"]["radesm"].as<float>();
    sphecamE1 =  config["Limiters"]["sphecamE1"].as<float>();
    sphecamI1 =  config["Limiters"]["sphecamI1"].as<float>();
    opamP1 =  config["Limiters"]["opamP1"].as<float>();
    opamR1 =  config["Limiters"]["opamR1"].as<float>();
    sphecamE2 =  config["Limiters"]["sphecamE2"].as<float>();
    sphecamI2 =  config["Limiters"]["sphecamI2"].as<float>();
    opamP2 =  config["Limiters"]["opamP2"].as<float>();
    opamR2 =  config["Limiters"]["opamR2"].as<float>();
    larsen_limiter = config["Limiters"]["larsen"].as<float>();
    leakArea1 = config["FixedParameters"]["leakArea1"].as<double>(); 
    leakArea2 = config["FixedParameters"]["leakArea2"].as<double>();
    leakMaterial = config["FixedParameters"]["leakMaterial"].as<int>();

    //RadIncomingFluxTemperature = config["BoundaryConditions"]["Rad_Flux_Temperature"].as<double>();
    RightFluidBoundaryCondition = config["BoundaryConditions"]["Hydro_Right_Boundary"].as<std::string>();
    RightRadBoundaryCondition = config["BoundaryConditions"]["Rad_Right_Boundary"].as<std::string>();
    LeftRadBoundaryCondition = config["BoundaryConditions"]["Rad_Left_Boundary"].as<std::string>();
    LeftDirichletEnergyDensity = config["BoundaryConditions"]["Rad_Dirichlet_Left"].as<double>();
    RightDirichletEnergyDensity = config["BoundaryConditions"]["Rad_Dirichlet_Right"].as<double>();
    LeftSourceTemperature =  config["BoundaryConditions"]["Rad_Source_Left"].as<double>();
    RightSourceTemperature =  config["BoundaryConditions"]["Rad_Source_Right"].as<double>();

    std::cout << "FIXED PARAMETERS INPUT " << "\n";
    std::cout << "NUMBER OF GRID CENTERS " << Nx <<"\n";
    std::cout << "NUMBER OF ENERGY GROUPS " << Ng <<"\n";
    std::cout << "RIGHT FLUID BOUNDARY CONDITION " << RightFluidBoundaryCondition << "\n";
    std::cout << "LEFT RAD BOUNDARY CONDITION " << LeftRadBoundaryCondition << "\n";
    std::cout << "RIGHT RAD BOUNDARY CONDITION " << RightRadBoundaryCondition << "\n";
    std::cout << "PRE HEAT STARTING INDEX " << PreHeatStartIndex <<"\n";
    std::cout << "PRE HEAT ENDING INDEX " << PreHeatLastIndex <<"\n";
    std::cout << "FRONT HEAT STARTING INDEX " << FrontHeatStartIndex <<"\n";
    std::cout << "FRONT HEAT ENDING INDEX " << FrontHeatLastIndex <<"\n";
    std::cout << "E HEAT FLUX LIMITER " << heatflxlmE <<"\n";
    std::cout << "I HEAT FLUX LIMITER " << heatflxlmI <<"\n";

    if (OutputFrequency == 0)
    {
        OutputFrequency = 1;
    }
    std::cout << std::endl;
}

// 