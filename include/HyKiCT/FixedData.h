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
#ifndef FIXEDDATA_HPP
#define FIXEDDATA_HPP
#include "Data.hpp"
#include "Switches.h"

class FixedData : public Data
{
    public:
        int Nx, Ng,RadNg, Np;
        double MaxE;
        float Cq = 2.0;
        double Gamma = 1.666667;
        int OutputFrequency;
        std::vector<double> Z, Ar, PreB, FrontB; 
        int PreHeatStartIndex = 0, PreHeatLastIndex = 0, FrontHeatStartIndex = 0, FrontHeatLastIndex = 0;
        int materialInterface;
        double ConstantOpacity;
        double leakArea1 = 0.0, leakArea2 = 0.0;
        int leakMaterial;
        //multipliers 
        float heatflxlmE = 0.0, heatflxlmI = 0.0, radflxlm = 1.000, prrm1 = 1.0, prrm2 = 1.0, radabm = 1.0, radbrm = 1.0, radcmm = 1.0, raddfm = 1.0, radesm = 1.0;
        float sphecamE1 = 1.0, sphecamI1 = 1.0,sphecamE2 = 1.0, sphecamI2 = 1.0, opamP1 = 1.0, opamR1 = 1.0,opamP2 = 1.0, opamR2 = 1.0;
        float radFFGauntB = 1.103;
        float larsen_limiter = 2.0;
        
        //Boundary Conditions
        std::string RightFluidBoundaryCondition = "f";
        std::string RightRadBoundaryCondition = "v";
        std::string LeftRadBoundaryCondition = "v";
        double RightDirichletEnergyDensity = 0 , LeftDirichletEnergyDensity = 0;
        double RightSourceTemperature, LeftSourceTemperature;
        std::vector<double> leakTimes,leakAreaMultiplier1,leakAreaMultiplier2;
        
        const double BOLTZMANN_CONSTANT = 1.38064852e-23;
        const double PROTON_MASS = 1.6726219e-27;
        const double ELECTRON_MASS = 9.10938356e-31;
        const double VACUUM_PERMITTIVITY = 8.85418782e-12;
        const double ELECTRON_CHARGE = 1.60217662e-19;
        const double HBAR = 1.054571800e-34;
        const double PLANCK_CONSTANT = 6.62607004e-34;
        const double STEFAN_BOLTZMANN_CONSTANT = 5.670374419e-08;
        const double SPEED_OF_LIGHT = 299792458;
        const double ev_to_k = ELECTRON_CHARGE / BOLTZMANN_CONSTANT;
        const double k_to_ev = BOLTZMANN_CONSTANT / ELECTRON_CHARGE;
        const double cm_to_m = 1E-2;
        const double m_to_cm = 1e2;
        const double me_mp_ratio = ELECTRON_MASS/PROTON_MASS;
        const double radiationConstant = 4 * STEFAN_BOLTZMANN_CONSTANT / SPEED_OF_LIGHT;
        

        void init(std::string path) override;
        void initVectors(std::string varFilePath, bool coupleMulti);
        void loadLeakAreas(std::string path);
};

#endif
