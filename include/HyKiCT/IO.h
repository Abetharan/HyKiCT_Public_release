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
#ifndef IO_HPP
#define IO_HPP
#include <map>
#include <yaml-cpp/yaml.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/stat.h>
class IO
{
    
    public:
        void createDirectories();
        std::vector<std::string> DumpVars;
        std::map<std::string, std::string> dir_paths;
        std::map<std::string, bool> dumpBools;
        std::vector<std::string> dumpVars{"VISCOSITY","VELOCITY","CELL_WALL_X","CELL_CENTRE_X","DENSITY","ION_NUMBER_DENSITY",
                                    "ELECTRON_NUMBER_DENSITY","ION_TEMPERATURE","ELECTRON_TEMPERATURE","ION_INTERNAL_ENERGY",
                                    "ELECTRON_INTERNAL_ENERGY","ELECTRON_PRESSURE","ION_PRESSURE","TOTAL_PRESSURE","ELECTRON_SPECIFIC_HEAT",
                                    "ELECTRON_DP_DT","ION_SPECIFIC_HEAT","ION_DP_DT","ELECTRON_HEAT_FLOW_X","ION_HEAT_FLOW_X",
                                    "ELECTRON_HEAT_CONDUCTION","ION_HEAT_CONDUCTION","INVERSE_BREM","EXCHANGE","PLASMA_FREQUENCY",
                                    "COULOMB_LOG", "TIME","LASER_STUFF","ZBAR", "RAD_FREE_FREE_EMISSION", "RAD_FREE_FREE_ABSORB",  "RAD_MULTI_GROUP",
                                    "RAD_TEMPERATURE", "RAD_ENERGY_DENSITY", "TOTAL_ION_SOURCE", "TOTAL_ELECTRON_SOURCE", "PDVWORK", "PLANCKOPACITIES", "ROSSLANDOPACITIES"};
        bool CreateOutputFolder;
        IO(){};
        IO(std::string path);    
        std::string FEOSPathMaterial1;
        std::string FEOSPathMaterial2;
        std::string OpacityPathMaterial1;
        std::string OpacityPathMaterial2;
        std::string InitPath;
        std::string OutPath;
        std::string SwitchPath;
        std::string LaserProfilePath;
        std::string LeakProfilePath;

};
#endif 