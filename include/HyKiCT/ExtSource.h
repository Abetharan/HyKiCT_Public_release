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
#ifndef EXTSOURCES_HPP
#define EXTSOURCES_HPP
#include <numeric>
#include <algorithm>
#include "GridData.h"
#include "FixedData.h"
#include "PlasmaParameters.h"
// #include "Switches.hpp"
class ExtSource
{
    private:
        std::vector<double> mLaserPowerOverTime, mTime;
        std::vector<double> mOpticalDepth, mCoulombLog; 
        double mLaserPower;
        double mLaserAngle = 0;
        bool mNoFullDump = false;
        double mAlbedo = 0;
        double mCriticalLeak = 0;
        double mCriticalDumpFraction = 1;
        int mCriticalSurface = -1;
        int mStartIndex = -1;
        std::vector<double> mBeta;
        void mCalculateCoulombLog(GridData *gridData);
        double mInterpolatePower(double time, double lower_time, double upper_time, double lower_power, double upper_power);
    public:
        // void InverseBrem(GridData *gridData, FixedData const &fixedData, Switches const &switchData);
        int mNx;
        double LaserDuration;
        double LaserWavelength;
        ExtSource(){};
        ExtSource(std::string path);
        void calculateBeta(GridData *gridData);
        void inverseBrem(GridData *gridData, FixedData const &fixedData);
        void invertedInverseBrem(GridData *gridData);
        void inverseBrem(GridData *gridData);
        void setNoFullDump(bool on);
        double getLaserPower();
        void updateLaserPower(double time); //For variable lasing
        void setLaserPower(double power);
        void setLaserWavelength(double wavelength);
        void setLaserPowerOverTime(std::vector<double> laserPower, std::vector<double> time);
        void Init(std::string path);
        void setBeta(std::vector<double> beta);
        void setStartIndex(int nx){mStartIndex = nx;};
};
#endif