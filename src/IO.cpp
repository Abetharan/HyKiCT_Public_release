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
#include "HyKiCT/IO.h"

IO::IO(std::string path)
{

    YAML::Node config = YAML::LoadFile(path);

    InitPath = config["Paths"]["Init_Path"].as<std::string>();
    OutPath = config["Paths"]["Out_Path"].as<std::string>();
    FEOSPathMaterial1 = config["Paths"]["FEOS_Path_1"].as<std::string>();
    FEOSPathMaterial2 = config["Paths"]["FEOS_Path_2"].as<std::string>();
    OpacityPathMaterial1 = config["Paths"]["Opacity_Path_1"].as<std::string>();
    OpacityPathMaterial2 = config["Paths"]["Opacity_Path_2"].as<std::string>();
    LaserProfilePath = config["Paths"]["Laser_Profile_Path"].as<std::string>();
    LeakProfilePath = config["Paths"]["Leak_Profile_Path"].as<std::string>();
    std::cout << "IO PARAMETERS INPUT" << "\n";
    std::cout << "INITIALISATION PATH: " << InitPath<<"\n";
    std::cout << "OUTPUT PATH: " << OutPath<<"\n";
    std::cout << "FEOS MATERIAL 1 PATH: " << FEOSPathMaterial1<<"\n";
    std::cout << "FEOS MATERIAL 2 PATH: " << FEOSPathMaterial2<<"\n";

    std::cout << "DUMPING FOLLOWING OUTPUTS" << "\n";
    std::cout << std::endl;
    bool coord = config["Output"]["Coordinates"].as<bool>();
    bool viscosity = config["Output"]["Viscosity"].as<bool>();
    bool vel = config["Output"]["Velocity"].as<bool>();
    bool heat = config["Output"]["HeatFlow"].as<bool>();
    bool ex = config["Output"]["Exchange"].as<bool>();
    bool ib = config["Output"]["InverseBrem"].as<bool>();
    bool ff = config["Output"]["RadFreeFree"].as<bool>();
    bool rt = config["Output"]["RadTemperature"].as<bool>();
    bool red = config["Output"]["RadEnergyDensity"].as<bool>();
    bool rmged = config["Output"]["RadMultiGroupEnergy"].as<bool>();
    bool p = config["Output"]["Pressure"].as<bool>();
    bool t = config["Output"]["Temperature"].as<bool>();
    bool ie = config["Output"]["InternalEnergy"].as<bool>();
    bool rho = config["Output"]["MassDensity"].as<bool>();
    bool n = config["Output"]["NumberDensity"].as<bool>();
    bool zbar = config["Output"]["Zbar"].as<bool>();
    bool pp = config["Output"]["PlasmaParameters"].as<bool>();
    CreateOutputFolder = config["Output"]["CreateOutputFolder"].as<bool>();
    bool time = true;
    
    std::vector<bool>bool_list{viscosity, vel, coord, coord, rho, n, n, t, t, ie, ie, p, p, p,
                                 pp, pp, pp ,pp, heat, heat, heat, heat, ib, ex, pp, pp, time, ib, zbar, ff,ff, rmged, rt, red, ex, ex,ex, rt,rt};
    for(unsigned long int i = 0; i < bool_list.size(); i++)
    {
        dumpBools[dumpVars[i]] = bool_list[i];
    }
}
    

void IO::createDirectories()
{
    std::string base_folder_path;
    for(unsigned long int i =0; i < dumpVars.size(); i++)
    {
        if(CreateOutputFolder) 
        {
            base_folder_path = OutPath + "/HyKiCT_OUTPUT/";
            if(i == 0)
            {
                mkdir(base_folder_path.c_str(), 0777);
            }
        }
        else 
        {
            base_folder_path = OutPath;
        }
        if(dumpBools[dumpVars[i]])
        {
            std::string dir_path = base_folder_path + '/' + dumpVars[i] + '/';
            std::cout << dir_path << '\n';
            dir_paths[dumpVars[i]] = dir_path;
            mkdir(dir_path.c_str(), 0777);
        }
        // if(status != 0)
        // {
        //     throw std::invalid_argument("PATH NOT CREATED");
        // }
    }
}