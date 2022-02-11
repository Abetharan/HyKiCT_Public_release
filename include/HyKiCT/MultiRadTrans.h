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
#ifndef MULTIRADTRANS_H
#define MULTIRADTRANS_H
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <numeric>
#include <optional>
#include <tuple>
#include "FixedData.h"
#include "GridData.h"
#include "Switches.h"
#include "TimeStep.h"
#include "LoadInTables.h"
// #include "MatrixSolver.h"
#include "RadTrans.hpp"
// #include "petscksp.h" 
// #include "petscmat.h"
// #include "petscvec.h"

class MultiRadTrans : public RadTrans
{
    private:
        int mNx;
        int mNg = 1;
        float mgff;
        double mLeakArea = 0;
        std::vector<double> mLeakAreaVector;
        int mMaterialInterface = 0;
        int mLeakMaterial = 0; 
        double mLocalEmissionCap = 1e-15;
        double mLocalEnergyDensityCap = 1e-25;
        double mLocalPlanckOpacityCap = 1e-25;
        double mLocalRossOpacityCap = 1e-25;
        //PetscInitialize();
        // PetscScalar mvalue[3],mLeftBoundaryValue[2],mRightBoundaryValue[2], mlastValue, mstartValue, mraddens;
        //Create PETSC vectors
        std::tuple<double,double> mLeftBoundaryValue{-1,-1};
        std::tuple<double, double> mRightBoundaryValue{-1,-1};
        std::vector<double> mEnergyDensities, mOpticalDepth;
        double mRightBoundaryEnergyDensity, mLeftBoundaryEnergyDensity;
        double mPlanckPiLFunction(double x, int l,bool without_1);
        double mReturnGamma(double x);
        void mLeftBoundaryCondition(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData);
        void mRightBoundaryCondition(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData);
        void mLeftReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d);
        void mRightReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d);
        void mImplicitLeftVacuumSourceBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in, double *a, double *b, double *c, double*d);
        void mLeftVacuumSourceBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in);
        void mRightVacuumSourceBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double f_in);
        void mSolveTridiagonal(double* a, double* b, double* c, double* d, int n);
        void mImplicitRightReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d);
        void mImplicitLeftReflectiveBoundaryCondition(GridData const *gridData, FixedData const&fixedData, TimeStep const &timeData, double *a, double *b, double *c, double*d,double heating );
        double mInterpolateArea(double time, double lower_time, double upper_time, double lower_power, double upper_power);
        void mFindLeakArea(double search_time, FixedData const &fixedData);
        // fill every step, and clear every step. 
        struct LocalData{
            std::vector<double> RadPlanckAbsorptionOpacity, 
                                RadPlanckEmissionOpacity, RadRossAbsorptionOpacity,
                                RadEnergyDensity, RadDiffusionCoefficient,RadFFEmission,RadFFAbsorb;
        };
        LocalData mLocal;
        int mEnergyGroup = 0;
        std::vector<double> mPlanck, mPlanckIntegral;
        std::vector<double> mLeakSourceCoefficient;
    public:
        std::vector<LoadInTables> mOpacityTables;
        double exp2(double x);
        double RightBoundaryValue, LeftBoundaryValue;
        MultiRadTrans(int nx);
        void setMaterialInterface(int index) override;
        void storeOpacityTables(std::vector<LoadInTables> table) override;
        void initEnergyDensity(GridData const *gridData, FixedData const &fixedData);
        void setPlanckianOpacity(double value);
        void setRossOpacity(double value);
        void setRadEnergyDensity(double value);
        void setRadEnergyDensity( std::vector<double> values);
        void setIntegratedRadEnergyDensity(std::vector<double> totalRadEnergyDensity);
        void setEnergyGroup(int group) override;
        void setTotalEnergyGroup(int totalGroups);
        void setPlanckIntegral(std::vector<double> &values);
        void setDiffusionCoefficient(std::vector<double> &values);
        void setInitLeakParams(int leakMaterial, double area);
        void setLeakArea(FixedData const &fixedData, double search_time);
        void getOpacities(GridData *gridData);
        void setAnalyticalOpacity(GridData const *gridData, FixedData const &fixedData);
        std::vector<double> getPlanckIntegral();
        std::vector<double> getLocalEnergyDensity();
        void heatRegion(std::vector<double> heating,TimeStep const &timeData);
        void fullyImplicit(GridData const* gridData, FixedData const& fixedData,TimeStep const &timeData);
        void fullyImplicit(GridData const* gridData, FixedData const& fixedData,TimeStep const &timeData, std::vector<double> heating);
        void LeakSource(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData);
        void calculateRadPressure(GridData *gridData, FixedData const &fixedData);
        // Free-Free Processes and its associated opacity calculations
        void calcFreeFreeAbsorption(GridData const *gridDataT,FixedData const&fixedData);
        void calcFreeFreeEmission(GridData const *gridDataT,FixedData const&fixedData);
        void planckIntegral(GridData const *gridDataT, FixedData const&fixedData);
        // Diffusion
        void ficksDiffusionCoefficient(GridData const *gridData, FixedData const &fixedData);
        void ficksDiffusion(GridData const *gridData, FixedData const &fixedData, TimeStep const &timeData);
        //Opacities
        void calculateOpacities(GridData const *gridData, FixedData const &fixedData,  bool skipSearch = false);
        //For Opacity Tables
        void planckianOpacity(GridData const *gridData,int const table, int const i, float const opamP);
        // void planckianEmissionOpacity(GridData *gridData,int table, int i);
        void rosslandOpacity(GridData const *gridData,int const table, int const i, float const opamR);
        //Radiation Energy Density update
        void updateRadiationEnergyDensity(GridData const *gridDataT,FixedData const&fixedData, TimeStep const &timeData);
        //misc
        void updateRadiationQuantities(GridData const *gridDataT, FixedData const &fixedData, TimeStep const &timeData, int updateType) override;
        void updateRadiationQuantities(GridData const *gridDataT, FixedData const &fixedData, TimeStep const &timeData);
        //non-constant
        void gather(GridData *gridData) override;
        void updateRadiationTemperature(GridData *gridData, FixedData const &fixedData);
        void copySearchEntries(MultiRadTrans &tables);
        std::vector<LoadInTables>& getSearchEntries();
};
#endif