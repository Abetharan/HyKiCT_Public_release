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
#ifndef SNBSOURCE_HPP
#define SNBSOURCE_HPP
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
// #include "petscksp.h" 
// #include "petscmat.h"
// #include "petscvec.h"
#include "GridData.h"
#include "FixedData.h"
#include "Switches.h"
#include "TimeStep.h"
#include "FastMaths.hpp"
// #include "MatrixSolver.h"

class SNBSource
{
    private:
        int mNx, mNg = 0, mGroupID;
        double mR = 1;
        bool mSeperated = true;
        std::tuple<double,double> mLeftBoundaryValue{-1,-1};
        std::tuple<double, double> mRightBoundaryValue{-1,-1};
        std::vector<double>  mDqNL, mEField,mSNBGridWall, mSNBGridCentered, mSpitzer;
        std::vector<std::vector<double>> mLambdaStarWall, mLambdaStarCentre, mLambdaE, 
                                        mSNBBetaWall, mSNBDbetaWall, mSNBBetaCentre, mSNBDbetaCentre,
                                        mSNBH, mSNBU;
        void mCalculateSNBBeta(std::vector<double> const &TemperatureE);
        void mCalculateSNBLambdaStar(std::vector<double> const &TemperatureE,
                                    std::vector<double> const &NumberDensityE, 
                                    std::vector<double> const &Zbar, 
                                    std::vector<double> const &CoulombLog, FixedData const &fixedData);
        void mCalculateSpitzerElectricField(std::vector<double> const &cell_centre_grid, 
                                            std::vector<double> const &TemperatureE, 
                                            std::vector<double> const &NumberDensityE, 
                                            std::vector<double> const &Zbar, FixedData const &fixedData);
        void mCalculateSNBLambdaE(FixedData const &fixedData);
        void mCalculateSNBU(std::vector<double> const &SpitzerHarmHeatFlow);
        // void mCalculateSNBH(std::vector<double> const &cell_centre_grid,
        //                      std::vector<double> const &cell_wall_grid,
                            //   std::vector<double> const &Zbar, MatrixSolver &matrixSolver);
        void mCalculateSNBH( GridData const *gridData,
                              std::vector<double> const &Zbar);
        void mCalculateSNBDqNL(std::vector<double> &SNBDqNL, std::vector<double> const &cell_centre_grid);
        void mSolveTridiagonal(double* a, double* b, double* c, double* d, int n);
    public:
        SNBSource(int nx);
        void InitSNB(int ng, int id, double maxE);
        void snbHeatFlowCorrection(GridData const *gridData, FixedData const &fixedData);
        void SNBGather(GridData *gridData);
        void SNBGatherH(GridData *gridData);
        void SNBCorrect(GridData *gridData);
        void switchToAveragedSNB(){mSeperated = false;};
};
#endif