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
#include "HyKiCT/ExtSource.h"
ExtSource::ExtSource(std::string path)
{
    YAML::Node config = YAML::LoadFile(path);
    mNx = config["FixedParameters"]["nx"].as<int>();
    mLaserPower = config["LaserParams"]["Constant_Power"].as<double>();
    LaserWavelength = config["LaserParams"]["Wavelength"].as<double>();
    mLaserAngle = config["LaserParams"]["Angle"].as<double>();
    LaserDuration = config["LaserParams"]["Duration"].as<double>();
    mAlbedo = config["LaserParams"]["Albedo"].as<double>();
    mCriticalLeak = config["LaserParams"]["CriticalLeak"].as<double>();
    mCriticalDumpFraction = config["LaserParams"]["CriticalDumpFraction"].as<double>();
    mStartIndex = config["LaserParams"]["StartIndex"].as<int>();
    if(mStartIndex < 0)
    {
        mStartIndex = mNx;
    }
    if(mCriticalDumpFraction < 1.0)
    {
        mNoFullDump = true;
    }
}
void ExtSource::Init(std::string path)
{
    std::cout << "Loading Laser Profile" << std::endl;
    std::ifstream initFile(path);
    std::string line;

    if (initFile.is_open())
    {
        double time;
        double intensity;
        
        while (getline(initFile, line))
        {
            std::istringstream ss(line);
            ss >> time >> intensity;
            mLaserPowerOverTime.push_back(intensity);
            mTime.push_back(time);
        }
    }
    LaserDuration = mTime[mTime.size() - 1];
}
void ExtSource::setLaserPower(double power)
{
    mLaserPower = power;
}
void ExtSource::setLaserWavelength(double wavelength)
{
    LaserWavelength = wavelength;
}
void ExtSource::setLaserPowerOverTime(std::vector<double> laserPower, std::vector<double> time)
{
    mLaserPowerOverTime = laserPower;
    mTime = time;
}
double ExtSource::mInterpolatePower(double time, double lower_time, double upper_time, double lower_power, double upper_power)
{
    double weight = (time - lower_time) / (upper_time - lower_time);
    double interIntensity = lower_power * (1 - weight) + upper_power * weight; 
    return interIntensity;
}
void ExtSource::updateLaserPower(double time)
{
    // int laserIndex = std::round(time / mTimeBin);
    int i = 0;
    double tol = time * 1e-5;
    bool noIndex = true; 
    while(noIndex)
    {
        if((time > mTime[mTime.size() - 1]) || (time < mTime[0]) || (time < 1e-30))
        {
            mLaserPower = 0;
            break;
        }
        else
        {
            auto time_upper_bound_search = std::lower_bound(mTime.begin(), mTime.end(), time );
            auto time_lower_bound_search = time_upper_bound_search - 1;
            int time_lower_index = time_lower_bound_search - mTime.begin();
            int time_upper_index = time_upper_bound_search - mTime.begin();
            mLaserPower = mInterpolatePower(time, mTime[time_lower_index], mTime[time_upper_index],mLaserPowerOverTime[time_lower_index], mLaserPowerOverTime[time_upper_index]);
            break;
            //Old Bin Method
            // auto time_iterator = std::lower_bound(mTime.begin(), mTime.end(), time);
            // int index = time_iterator - mTime.begin();
            // if(time < mTime[index])
            // {
            //     index--;
            // }
            // mLaserPower = mLaserPowerOverTime[index]*std::cos(M_PI*mLaserAngle/180);
            // break;
        }
    }
}
double ExtSource::getLaserPower()
{
    return mLaserPower;
}
void ExtSource::setNoFullDump(bool on)
{
    mNoFullDump = on;
}

//REF LEGACY Remove at some point
// void ExtSource::inverseBrem(GridData *gridData, FixedData const &fixedData)
// {
//     int nx = fixedData.Nx;
//     double criticalDensity = 1114326918632954.5 / std::pow(LaserWavelength, 2);
//     mOpticalDepth.resize(nx+1, 0); 
//     double dyinglaser = mLaserPower;  
//     int critical_surface = -1;
//     double DistanceTravelled;
//     double alpha;
//     for(int i = nx; i >= 0; i--)
//     {
//         if(i == nx)
//         {
//             continue;
//         }
//         DistanceTravelled = abs(gridData->CellWallCoord[i+1] - gridData->CellWallCoord[i]);
//         if(i == 0)
//         {
//             double beta = gridData->InterNumberDensityE[0] /criticalDensity; 
//             if(beta >= 1)
//             {
//                 critical_surface = i;
//                 break;
//             }
//             alpha = 1.2130370535544696E-14 * gridData->InterNumberDensityE[0] * gridData->InterZbar[0] * gridData->InterCoulombLog[0] 
//                             * std::pow(gridData->InterTemperatureE[0], -1.5) * std::pow(beta, 2) * std::pow((1 - beta), -0.5);
//         }
//         else
//         {
//             double beta = gridData->InterNumberDensityE[2*i - 1] /criticalDensity; 
//             if(beta >= 1)
//             {
//                 critical_surface = i;
//                 break;
//             }
//             alpha = 1.2130370535544696E-14 * gridData->InterNumberDensityE[2*i - 1] * gridData->InterZbar[2*i - 1] * gridData->InterCoulombLog[2*i -1] 
//                             * std::pow(gridData->InterTemperatureE[2*i - 1], -1.5) * std::pow(beta, 2) * std::pow((1 - beta), -0.5);
//         }

//         mOpticalDepth[i] = alpha * DistanceTravelled;       
//     }
//     double totalPower = 0;
//     for(int i = nx - 1; i >=0; i--)
//     {
//         if(i == critical_surface)
//         {
//             gridData->TransmittedLaser[i] = 0;
//             double accumulatedPower = std::accumulate(gridData->PowerAbsorbed.begin(), gridData->PowerAbsorbed.end(), 0.0);
//             gridData->PowerAbsorbed[i] = mLaserPower - accumulatedPower;
//             gridData->InverseBrem[i] = (1/gridData->Mass[i]) * gridData->PowerAbsorbed[i];
//             break;
//         }
//         else
//         {
//             gridData->TransmittedLaser[i] = mLaserPower * std::exp(-1*mOpticalDepth[i]);
//             gridData->PowerAbsorbed[i] = mLaserPower * abs(std::exp(-1 * mOpticalDepth[i+1]) -
//                                                     std::exp(-1 * mOpticalDepth[i]));
//             // dyinglaser = gridData->TransmittedLaser[i];
//             if(gridData->PowerAbsorbed[i] < 0)
//             {
//                 std::cout << "Negative Energy " << std::endl;
//             }
//             totalPower+=gridData->PowerAbsorbed[i];
//             if (totalPower > mLaserPower)
//             {
//                 double diff_max = abs(mLaserPower - totalPower); 
//                 gridData->PowerAbsorbed[i] -=  diff_max;
//                 gridData->InverseBrem[i] = (1/gridData->Mass[i]) * gridData->PowerAbsorbed[i];
//                 break;
//             }
//         }
//         gridData->InverseBrem[i] = (1/gridData->Mass[i]) * gridData->PowerAbsorbed[i];
//     }
//     gridData->TransmittedLaser[nx] = mLaserPower;
//     double accumulatedPower = std::accumulate(gridData->PowerAbsorbed.begin(), gridData->PowerAbsorbed.end(), 0.0);
//     double relativePower = abs(accumulatedPower - mLaserPower) / mLaserPower;
//     if((!mNoFullDump) && (mLaserPower > 0))
//     {
//         assert(relativePower < 1e-3);
//     }
// }


void ExtSource::calculateBeta(GridData *gridData)
{
    int nx = mNx;
    double criticalDensity = 1114326918632954.5 / std::pow(LaserWavelength, 2);
    mBeta.resize(nx + 1, 0);
    for(int i = mStartIndex; i >= 0; i--)
    {
        if(i == 0)
        {
            mBeta[i] = 0.0;
        }
        else if(i == mNx)
        {
            // double beta = gridData->NumberDensityE[mNx - 1]/ criticalDensity;
            mBeta[i] = 0.0;//pow((beta), 1) * pow(1 - (beta), -0.5);
        }
        else 
        {
            double cell_wall_ne1 = (gridData->NumberDensityE[i] + gridData->NumberDensityE[i-1]) / 2;
            // double cell_wall_ne2 = (gridData->NumberDensityE[i] + gridData->NumberDensityE[i+1]) / 2;
            double cell_wall_ne = 2 * (gridData->NumberDensityE[i] * gridData->NumberDensityE[i - 1]) / (gridData->NumberDensityE[i] + gridData->NumberDensityE[i - 1]);
            double walled_beta = cell_wall_ne / criticalDensity ; 
            if (walled_beta > 1)
            {
                mCriticalSurface = i;
                break;
            }
            mBeta[i] = pow((walled_beta), 1) * pow(1 - (walled_beta), -0.5); 
            
        }
    }
}
void ExtSource::setBeta(std::vector<double> beta)
{
    mBeta = beta;
}
void ExtSource::inverseBrem(GridData *gridData)
{
    double DistanceTravelled, alpha;
    for(int i = mStartIndex; i>=0; i--)
    {
        if(i == mCriticalSurface)
        {
            break;
        }
        if(i == mStartIndex)
        {
            // intensity[i] = mLaserPower * (1 - mAlbedo)
            gridData->TransmittedLaser[i] = mLaserPower * (1 - mAlbedo);
        }
        else
        {   
            DistanceTravelled = abs(gridData->CellWallCoord[i+1] - gridData->CellWallCoord[i]); 
            //Linear 
            // alpha = 1.2130370535544696E-14 * gridData->NumberDensityE[i] * gridData->Zbar[i] * gridData->CoulombLogEI[i] 
            //                 * std::pow(gridData->TemperatureE[i], -1.5) * ((mBeta[i + 1] + mBeta[i])/2);
            alpha = gridData->CollisionFrequencyEIOverC[i] * ((mBeta[i + 1] + mBeta[i])/2);
            //Quadratic
            // double component = (beta[i] + 4*(pow((gridData->NumberDensityE[i] / criticalDensity), 2) * pow(1 - (gridData->NumberDensityE[i] / criticalDensity), -0.5)) +beta[i+1]) * DistanceTravelled/6;
            // alpha = 1.2130370535544696E-14 * gridData->NumberDensityE[i] * gridData->Zbar[i] * gridData->CoulombLogEI[i] 
                            // * std::pow(gridData->TemperatureE[i], -1.5) * component;
            // intensity[i] = intensity[i+1] * exp(-1 * alpha * DistanceTravelled);
            gridData->TransmittedLaser[i] = gridData->TransmittedLaser[i + 1] * exp(-1 * alpha * DistanceTravelled);//intensity[i];
        }
    }
    std::vector<double> diff(mNx, 0);
    diff = gridData->TransmittedLaser; 
    // std::adjacent_difference(gridData->TransmittedLaser.begin(), gridData->TransmittedLaser.end(),gridData->TransmittedLaser.begin());
    for(int i = mStartIndex; i>0; i--)
    {
        gridData->PowerAbsorbed[i - 1] = gridData->TransmittedLaser[i] - gridData->TransmittedLaser[i - 1];//intensity[i] - intensity[i - 1];
        
        if(gridData->PowerAbsorbed[i - 1] < 0)
        {
            std::cout << "Negative Energy " << std::endl;
        }
        if(gridData->PowerAbsorbed[i - 1] < 1e-30)
        {
            gridData->PowerAbsorbed[i - 1] = 0.0;
        }
        gridData->InverseBrem[i - 1] = gridData->PowerAbsorbed[i - 1] / gridData->Mass[i - 1];
    }
    //Do some modifications to critical 
    if((mCriticalDumpFraction < 1.0) && (mCriticalSurface>0))
    {
        gridData->PowerAbsorbed[mCriticalSurface] *= mCriticalDumpFraction;
        gridData->InverseBrem[mCriticalSurface] = gridData->PowerAbsorbed[mCriticalSurface] / gridData->Mass[mCriticalSurface];
    }
    if((mCriticalLeak > 0) && (mCriticalSurface>0))
    {
        double leak = gridData->PowerAbsorbed[mCriticalSurface]*mCriticalLeak;
        gridData->PowerAbsorbed[mCriticalSurface] -= leak;
        gridData->PowerAbsorbed[mCriticalSurface - 1] = leak;
        gridData->InverseBrem[mCriticalSurface] = gridData->PowerAbsorbed[mCriticalSurface] / gridData->Mass[mCriticalSurface];
        gridData->InverseBrem[mCriticalSurface - 1] = gridData->PowerAbsorbed[mCriticalSurface - 1] / gridData->Mass[mCriticalSurface - 1];
    }

    double accumulatedPower = std::accumulate(gridData->PowerAbsorbed.begin(), gridData->PowerAbsorbed.end(), 0.0);
    double relativePower = abs(accumulatedPower - mLaserPower*(1 - mAlbedo)) / (mLaserPower*(1 - mAlbedo));    
    if((!mNoFullDump) && (mLaserPower > 0))
    {
        assert(relativePower < 1e-3);
    }
    mCriticalSurface = -1;
}
