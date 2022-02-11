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
#include <algorithm>
#include <vector>
#include <utility>
#include "HyKiCT/GridData.h"
#include "HyKiCT/FixedData.h"
#include "HyKiCT/FluidDynamics.h"
// #include "HyKiCT/GrayRadTrans.h"
#include "HyKiCT/MultiRadTrans.h"
#include "HyKiCT/RadTrans.hpp"
#include "HyKiCT/TimeStep.h"
#include "HyKiCT/IO.h"
#include "HyKiCT/IdealGasEoS.h"
#include "HyKiCT/IntSource.h"
// #include "HyKiCT/MatrixSolver.h"
#include "HyKiCT/Init.h"
#include "HyKiCT/LoadInTables.h"
// #include <petscksp.h>

class ImplicitGrayRadiationTransportEquil:public::testing::TestWithParam<std::tuple<double>>
{
};
// void InitRadSources(GridData *gridData, double T_0, double R, int nx);
// void InitRadDiffusion(GridData *gridData, double U, double D_0, double R, int nx);