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
#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "GridData.h"
#include "FixedData.h"
#include "Data.hpp"
#include "IO.h"
class TimeStep : Data
{
    private:
        float mGamma = 1.66666666666667;
        int mNx, mAdaptiveMethod = 0;
        double stabilityConstraint(std::vector<double> lastVar, std::vector<double> currVar);
        double courantStabilityConstraint(std::vector<double> currDen,std::vector<double> currPressure, std::vector<double> currCoord);
        
        double mPVLimiter = 0.01, mVolumeLimiter = 0.01, mRadiationLimiter = 0.005, mElectronTemperatureLimiter = 0.005, mIonTemperatureLimiter =0.005, mCourantCondition = 0.85; 
        double mCalcSoundSpeed(double pressure, double density);
        double mCourant(double sound_speed, double dx);
        double mVolumePressureTimeLimit(double volume, double pressure, double dx);
        double mVolumeTimeLimit(double volumeT, double volumeT_1);
        double mRadiationTimeLimit(double RadiationT, double RadiationT_1);
        double mElectronTemperatureTimeLimit(double electron_temperartureT, double electron_temperartureT_1);
        double mIonTemperatureTimeLimit(double ion_temperartureT, double ion_temperartureT_1);

    public:
        TimeStep(int nx, float gamma);
        TimeStep(int nx);
        
        double Dt1, Dt05, InitialDt, dtGlobalMax, dtGlobalMin;
        double TotalTime;
        double Tmax;
        int MaxSteps;
        int Step = 0;
        void setTimeParameters(double pvlimit, double vollimit, double radlimit, double telimit,double tilimit, double cfl); 
        void setDt1(double dt1);
        void setDt05(double dt05);
        void setTotalTime();
        void calculateNewTime(GridData *gridDataT_1,  GridData *gridDataT);
        void incrementStep(){Step++;};
        std::string fileNameCreator(std::string path, std::string var, std::string step, std::string extension);

        void init(std::string path) override;
        void writeOut(std::ofstream &variable, std::vector<double> array, int length) override;
        void dump(IO &ioData) override;
};
#endif