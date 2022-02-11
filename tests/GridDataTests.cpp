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
#include <cstdlib>

class GridDataInitTests:public ::testing::TestWithParam<std::tuple<int>>
{
    protected:
        GridData *gridData;
        std::string loadPath = getenv("TEST_PATH");
        void TearDown()
        {
            delete gridData;
        }
};

TEST_P(GridDataInitTests, LoadCheck)
{
    int nx = std::get<0>(GetParam());
    std::string nx_s =std::to_string(nx);
    std::string localPath = loadPath + "/LOAD_IN_HEAT_COMPONENT_CHECK/" + nx_s;
    std::cout << localPath << '\n';
    gridData = new GridData(nx, true);
    gridData->init(localPath);
    EXPECT_EQ(gridData->Velocity.size(), nx + 1);
    EXPECT_EQ(gridData->CellWallCoord.size(), nx + 1);
    EXPECT_EQ(gridData->TemperatureE.size(), nx);
    EXPECT_EQ(gridData->TemperatureI.size(), nx);
    EXPECT_EQ(gridData->Density.size(), nx);
}

TEST_P(GridDataInitTests, LoadCheckMultiCheck)
{
    int nx = std::get<0>(GetParam());
    std::string nx_s =std::to_string(nx);
    std::string localPath = loadPath + "/LOAD_IN_HEAT_COMPONENT_CHECK/" + nx_s;
    gridData = new GridData(nx, true);
    gridData->init(localPath, false, true, 1);

    EXPECT_EQ(gridData->Velocity.size(), nx + 1);
    EXPECT_EQ(gridData->CellWallCoord.size(), nx + 1);
    EXPECT_EQ(gridData->TemperatureE.size(), nx);
    EXPECT_EQ(gridData->TemperatureI.size(), nx);
    EXPECT_EQ(gridData->Density.size(), nx);
    EXPECT_EQ(gridData->HeatFlowMultiplier.size(), nx);
}

TEST_P(GridDataInitTests, LoadCheckDivQCheck)
{
    int nx = std::get<0>(GetParam());
    std::string nx_s =std::to_string(nx);
    std::string localPath = loadPath + "/LOAD_IN_HEAT_COMPONENT_CHECK/" + nx_s;
    gridData = new GridData(nx, true);
    gridData->init(localPath, false, true, 2);
    EXPECT_EQ(gridData->Velocity.size(), nx + 1);
    EXPECT_EQ(gridData->CellWallCoord.size(), nx + 1);
    EXPECT_EQ(gridData->TemperatureE.size(), nx);
    EXPECT_EQ(gridData->TemperatureI.size(), nx);
    EXPECT_EQ(gridData->Density.size(), nx);
    EXPECT_EQ(gridData->LoadInDivQ.size(), nx);
}

INSTANTIATE_TEST_SUITE_P(
    CompleteInitCheck,
    GridDataInitTests,
    ::testing::Values(
        std::make_tuple(1E2),
        std::make_tuple(555),
        std::make_tuple(1E3),
        std::make_tuple(1E4)));