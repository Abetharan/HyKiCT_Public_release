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
#ifndef DATA_HPP
#define DATA_HPP
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <stdexcept>
#include <map>
#include "IO.h"
class Data
{
    public:
        virtual void dump(IO &ioData) {};
        virtual void init(std::string path) = 0;
        virtual void writeOut(std::ofstream &variable, std::vector<double> array, int length) {};
        virtual void terminalPrint() {};
};
#endif