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
#ifndef GRAYRADTRANS_H
#define GRAYRADTRANS_H
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <numeric>
#include <tuple>
#include "FixedData.h"
#include "GridData.h"
#include "Switches.h"
#include "TimeStep.h"
#include "LoadInTables.h"
#include "MatrixSolver.h"
#include "RadTrans.hpp"
#include "petscksp.h" 
#include "petscmat.h"
#include "petscvec.h"
class GrayRadTrans : public RadTrans
{
    private:
        int mNx;
        float mgff;
        int mMaterialInterface;
        std::tuple<double,double> mLeftBoundaryValue{-1,-1};
        std::tuple<double, double> mRightBoundaryValue{-1,-1};
        std::vector<double> mEnergyDensities, mOpticalDepth;
        std::vector<LoadInTables> mOpacityTables;
        double mRightBoundaryEnergyDensity, mLeftBoundaryEnergyDensity; 
        void mLeftBoundaryCondition(GridData *gridData, FixedData const &fixedData, TimeStep const &timeData);
        void mRightBoundaryCondition(GridData *gridData, FixedData const &fixedData, TimeStep const &timeData);
        void mLeftReflectiveBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData);
        void mRightReflectiveBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData);
        void mLeftVacuumSourceBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in);
        void mRightVacuumSourceBoundaryCondition(GridData *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in);

    public:
        double RightBoundaryValue, LeftBoundaryValue;
        GrayRadTrans(int nx);
        void setMaterialInterface(int index) override;
        void storeOpacityTables(std::vector<LoadInTables> table) override;
        void initEnergyDensity(GridData *gridData, FixedData const &fixedData);
        void setPlanckianOpacity(GridData *gridData, double value);
        void setRossOpacity(GridData *gridData, double value);
        void setRadEnergyDensity(GridData *gridData, double value);
        void setRadEnergyDensity(GridData *gridData, std::vector<double> values);
        
        // Radiation Pressure Effects
        void radationPressure(GridData *gridData);
        
        // Free-Free Processes and its associated opacity calculations
        void calcFreeFreeAbsorption(GridData *gridDataT,FixedData const&fixedData);
        void calcFreeFreeEmission(GridData *gridDataT,FixedData const&fixedData);
        
        // Diffusion
        void ficksDiffusionCoefficient(GridData *gridData, FixedData const &fixedData);
        void ficksDiffusion(GridData *gridData, FixedData const &fixedData, TimeStep const &timeData, MatrixSolver &matrixSolver);
        
        //Opacities
        void calculateOpacities(GridData *gridData, bool calculateRoss, bool skipSearch = false);
        //For Opacity Tables
        void planckianOpacity(GridData *gridData,int table, int i);
        // void planckianEmissionOpacity(GridData *gridData,int table, int i);
        void rosslandOpacity(GridData *gridData,int table, int i);

        //Radiation Energy Density update
        void updateRadiationEnergyDensity(GridData *gridDataT,FixedData const&fixedData, TimeStep const &timeData);
        void updateRadiationEnergyDensity(GridData *gridDataT);
        
        //misc
        void updateRadiationTemperature(GridData *gridData, FixedData const &fixedData);
        void updateRadiationQuantities(GridData *gridDataT, FixedData const &fixedData, TimeStep const &timeData, MatrixSolver &matrixSolver, int updateType) override;
};
#endif
