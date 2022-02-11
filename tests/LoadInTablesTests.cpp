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
#include "HyKiCT/GridData.h"
#include "HyKiCT/LoadInTables.h"
#include <iostream> 
#include <numeric>
#include <algorithm>
#include <array>
#include <string>

TEST(LoadTest, EoS)
{
    LoadInTables material_1(0, 1);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LoadInTablesTest/EoSTest";
    material_1.feosEoSLoadIn(localPath);
    // material_1.converTemperatureToeV(fixedData.ev_to_k);
    EXPECT_EQ (5, material_1.FEOSTemperature.size());
    EXPECT_EQ (6, material_1.FEOSDensity.size());
    EXPECT_EQ (30, material_1.intEnergyElectronFeosTable.size());
    EXPECT_EQ (30, material_1.pressureElectronFeosTable.size());
    EXPECT_EQ (30, material_1.intEnergyIonFeosTable.size());
    EXPECT_EQ (30, material_1.pressureIonFeosTable.size());
    ASSERT_EQ(30, material_1.pressureIonFeosTable.size());
}
TEST(LoadTest, Opacity)
{
    LoadInTables material_1(0, 1);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LoadInTablesTest/";
    material_1.OpacityLoadIn(localPath);
    EXPECT_EQ (5, material_1.OpacityTemperature.size());
    EXPECT_EQ (6, material_1.OpacityDensity.size());
    EXPECT_EQ (30, material_1.rosslandOpacity.size());
    EXPECT_EQ (30, material_1.planckEmissionOpacity.size());
}
TEST(searchTest, EoS)
{
    LoadInTables material_1(0, 1);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LoadInTablesTest/EoSTest";
    material_1.feosEoSLoadIn(localPath);
    // material_1.converTemperatureToeV(fixedData.ev_to_k);
    GridData initGridData(1, true);
    initGridData.TemperatureE = {15.0};
    initGridData.TemperatureI = {15.0};
    initGridData.Density = {0.1};
    material_1.findAllFeosIndicies(&initGridData);
    int elower_temp_index = std::get<0>(material_1.searchQuantsElectron[0].te_index);
    int eupper_temp_index = std::get<1>(material_1.searchQuantsElectron[0].te_index);
    int elower_density_index = std::get<0>(material_1.searchQuantsElectron[0].rho_index);
    int eupper_density_index = std::get<1>(material_1.searchQuantsElectron[0].rho_index);
    

    std::array<int, 4> trueCombinedindex{1,2,6,7};
    std::array<int, 4> trueTeSearchindex{1,2};
    std::array<int, 4> trueRhoSearchindex{0,1};
    EXPECT_EQ(trueCombinedindex[0], material_1.searchQuantsElectron[0].combined_index[0]);
    EXPECT_EQ(trueCombinedindex[1], material_1.searchQuantsElectron[0].combined_index[1]);
    EXPECT_EQ(trueCombinedindex[2], material_1.searchQuantsElectron[0].combined_index[2]);
    EXPECT_EQ(trueCombinedindex[3], material_1.searchQuantsElectron[0].combined_index[3]);
    EXPECT_EQ(trueTeSearchindex[0], elower_temp_index);
    EXPECT_EQ(trueTeSearchindex[1], eupper_temp_index);
    EXPECT_EQ(trueRhoSearchindex[0], elower_density_index);
    EXPECT_EQ(trueRhoSearchindex[1], eupper_density_index);
}
TEST(searchTest, Opacity)
{
    LoadInTables material_1(0, 1);
    std::string loadPath = getenv("TEST_PATH");
    std::string localPath = loadPath + "/LoadInTablesTest/";
    material_1.OpacityLoadIn(localPath);
    // material_1.converTemperatureToeV(fixedData.ev_to_k);
    GridData initGridData(1, true);
    initGridData.TemperatureE = {15.0};
    initGridData.TemperatureI = {15.0};
    initGridData.Density = {0.1};
    material_1.findAllOPACITYIndicies(&initGridData);
    int elower_temp_index = std::get<0>(material_1.searchQuantsOpacity[0].te_index);
    int eupper_temp_index = std::get<1>(material_1.searchQuantsOpacity[0].te_index);
    int elower_density_index = std::get<0>(material_1.searchQuantsOpacity[0].rho_index);
    int eupper_density_index = std::get<1>(material_1.searchQuantsOpacity[0].rho_index);
    

    std::array<int, 4> trueCombinedindex{1,2,6,7};
    std::array<int, 4> trueTeSearchindex{1,2};
    std::array<int, 4> trueRhoSearchindex{0,1};
    EXPECT_EQ(trueCombinedindex[0], material_1.searchQuantsOpacity[0].combined_index[0]);
    EXPECT_EQ(trueCombinedindex[1], material_1.searchQuantsOpacity[0].combined_index[1]);
    EXPECT_EQ(trueCombinedindex[2], material_1.searchQuantsOpacity[0].combined_index[2]);
    EXPECT_EQ(trueCombinedindex[3], material_1.searchQuantsOpacity[0].combined_index[3]);
    EXPECT_EQ(trueTeSearchindex[0], elower_temp_index);
    EXPECT_EQ(trueTeSearchindex[1], eupper_temp_index);
    EXPECT_EQ(trueRhoSearchindex[0], elower_density_index);
    EXPECT_EQ(trueRhoSearchindex[1], eupper_density_index);
}