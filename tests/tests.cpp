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
// #include <petscksp.h>
#include "RadiationTests.cpp"
#include "ImplicitRadiationTests.cpp"
#include "IntSourceTests.cpp"
#include "GridDataTests.cpp"
#include "ExtSourceTests.cpp"
#include "LoadInTablesTests.cpp"
#include "FluidDynamicsTests.cpp"
int main(int argc, char **argv)
{
    std::ostringstream local;
    auto cout_buff = std::cout.rdbuf(); // save pointer to std::cout buffer
    std::cout.rdbuf(local.rdbuf()); 
    ::testing::InitGoogleTest(&argc, argv);
    int re = RUN_ALL_TESTS();
    std::cout.rdbuf(cout_buff);
    return re;
}