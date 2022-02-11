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
#include "HyKiCT/TimeStep.h"


TimeStep::TimeStep(int nx, float gamma)
{
    mNx = nx;
    mGamma = gamma;
    TotalTime = 0.0;
}

TimeStep::TimeStep(int nx)
{
    mNx = nx;
    TotalTime = 0.0;
}
void TimeStep::setDt1(double dt1)
{
    if(dt1 < 0)
    {
        assert("Dt01 HAS TO BE GREATER THAN 0");
    }
    Dt1 = dt1;
}

void TimeStep::setDt05(double dt05)
{
    if(dt05 < 0)
    {
        assert("Dt05 HAS TO BE GREATER THAN 0");
    }
    Dt05 = dt05;
}

void TimeStep::setTotalTime()
{
    TotalTime += Dt05;
}

// void TimeStep::calculateNewTime(GridData *gridDataT_1,  GridData *gridDataT)
// {
    // //Based on Bucky, though courant conditions remained as traditional
    // std::vector<double> courant, pv, vol, rad, Te, Ti;
    // std::vector<double> time_conditions(dtGlobalMax);
    // time_conditions.reserve(7);
    // for(int i = 0; i < mNx; i++)
    // {
    //     double cs = mCalcSoundSpeed(gridDataT->TotalPressure[i], gridDataT->Density[i]);
    //     courant.push_back(mCourant(cs, gridDataT->Dx_k_1_2[i]));
    //     pv.push_back(mVolumePressureTimeLimit(1/gridDataT->Density[i], gridDataT->TotalPressure[i], gridDataT->Dx_k_1_2[i]));
    //     vol.push_back(mVolumeTimeLimit(1/gridDataT->Density[i], 1/gridDataT_1->Density[i]));
    //     rad.push_back(mRadiationTimeLimit(gridDataT->RadEnergyDensity[i], gridDataT_1->RadEnergyDensity[i]));
    //     Te.push_back(mElectronTemperatureTimeLimit(gridDataT->TemperatureE[i], gridDataT_1->TemperatureE[i]));
    //     Ti.push_back(mIonTemperatureTimeLimit(gridDataT->TemperatureI[i], gridDataT_1->TemperatureI[i]));
    // }
    // time_conditions.push_back(mCFL* *std::max_element(courant.begin(), courant.end()));
    // time_conditions.push_back(mPVLimiter* *std::max_element(pv.begin(), pv.end()));
    // time_conditions.push_back((mVolumeLimiter*Dt05) / *std::max_element(vol.begin(), vol.end()));
    // time_conditions.push_back((mRadiationLimiter*Dt05) / *std::max_element(rad.begin(), rad.end()));
    // time_conditions.push_back((mElectronTemperatureLimiter*Dt05) / *std::max_element(Te.begin(), Te.end()));
    // time_conditions.push_back((mIonTemperatureLimiter*Dt05) / *std::max_element(Ti.begin(), Ti.end()));

//     double min_time_conditions = *std::min_element(time_conditions.begin(), time_conditions.end());
//     if((min_time_conditions == INFINITY) || (min_time_conditions == 0))
//     {
//         min_time_conditions = Dt05;
//     }
//     std::vector<double> dt_vector = {dtGlobalMin, min_time_conditions};
//     double dt_3_2 = *std::max_element(dt_vector.begin(), dt_vector.end());
    
//     if((Tmax - TotalTime < dt_3_2) && Tmax > 0)
//     {
//         dt_3_2 = Tmax - TotalTime;
//     }
//     Dt1 = 0.5 * (dt_3_2 + Dt05);
//     Dt05 = dt_3_2;
// }
void TimeStep::calculateNewTime(GridData *gridDataT_1,  GridData *gridDataT)
{
    //Based on Bucky, though courant conditions remained as traditional
    std::vector<double> courant(mNx), pv(mNx), vol(mNx), rad(mNx), Te(mNx), Ti(mNx);
    std::vector<double> time_conditions(dtGlobalMax);
    time_conditions.reserve(7);
    for(int i = 0; i < mNx; i++)
    {
        double cs = mCalcSoundSpeed(gridDataT->TotalPressure[i], gridDataT->Density[i]);
        courant[i] = (mCourant(cs, gridDataT->Dx_k_1_2[i]));
        pv[i] = (mVolumePressureTimeLimit(1/gridDataT->Density[i], gridDataT->TotalPressure[i], gridDataT->Dx_k_1_2[i]));
        vol[i] = (mVolumeTimeLimit(1/gridDataT->Density[i], 1/gridDataT_1->Density[i]));
        rad[i] = (mRadiationTimeLimit(gridDataT->RadEnergyDensity[i], gridDataT_1->RadEnergyDensity[i]));
        Te[i] = (mElectronTemperatureTimeLimit(gridDataT->TemperatureE[i], gridDataT_1->TemperatureE[i]));
        Ti[i] = (mIonTemperatureTimeLimit(gridDataT->TemperatureI[i], gridDataT_1->TemperatureI[i]));
    }
    time_conditions.push_back(mCourantCondition* *std::max_element(courant.begin(), courant.end()));
    time_conditions.push_back(mPVLimiter* *std::max_element(pv.begin(), pv.end()));
    time_conditions.push_back((mVolumeLimiter*Dt05) / *std::max_element(vol.begin(), vol.end()));
    time_conditions.push_back((mIonTemperatureLimiter*Dt05) / *std::max_element(Ti.begin(), Ti.end()));
    if(mAdaptiveMethod == 0)
    {
        time_conditions.push_back((mRadiationLimiter*Dt05) / *std::max_element(rad.begin(), rad.end()));
        time_conditions.push_back((mElectronTemperatureLimiter*Dt05) / *std::max_element(Te.begin(), Te.end()));
    }
    else
    {
        std::vector<double> multiply_te, multiply_er;
        std::transform(Te.begin() + 1, Te.end(), 
                        gridDataT->TemperatureE.begin()+1, std::back_inserter(multiply_te),
                        std::multiplies<double>());
        std::transform(rad.begin() + 1, rad.end(), 
                        gridDataT->RadEnergyDensity.begin()+1, std::back_inserter(multiply_er),
                        std::multiplies<double>());
        
        double acc_multi_te = std::accumulate(multiply_te.begin(), multiply_te.end(), 0); 
        double acc_te = std::accumulate(gridDataT->TemperatureE.begin(), gridDataT->TemperatureE.end(), 0);
        double acc_multi_er = std::accumulate(multiply_er.begin(), multiply_er.end(), 0);
        double acc_er = std::accumulate(gridDataT->RadEnergyDensity.begin(), gridDataT->RadEnergyDensity.end(), 0);
        
        double dte = pow(acc_multi_te / acc_te, 0.5);
        double dter = pow(acc_multi_er / acc_er, 0.5);
        
        time_conditions.push_back(Dt05/dte * mElectronTemperatureLimiter);  // usually mElectronTemperatureLimiter is .1
        time_conditions.push_back(Dt05/dter * mRadiationLimiter); // usually mRadiationLimiter has vlaue 0.0005 

    }
    for(auto &i:time_conditions) //This accounts for any case where all the relative changse are negative
    {
        if(i < 0)
        {
            i*= -1; 
        }
        if(i == 0)
        {
            i = 1e30;
        }
    }
    double min_time_conditions = *std::min_element(time_conditions.begin(), time_conditions.end());
    if((min_time_conditions == INFINITY) || (min_time_conditions == 0))
    {
        min_time_conditions = Dt05;
    }
    std::vector<double> dt_vector = {dtGlobalMin, min_time_conditions};
    double dt_3_2 = *std::max_element(dt_vector.begin(), dt_vector.end());

    if((Tmax - TotalTime < dt_3_2) && Tmax > 0)
    {
        dt_3_2 = Tmax - TotalTime;
    }
    if(dt_3_2> Dt05*1.2)
    {
        dt_3_2 = Dt05*1.2;
    }
    if(dt_3_2 > dtGlobalMax)
    {
        dt_3_2 = dtGlobalMax;
    }
    Dt1 = 0.5 * (dt_3_2 + Dt05);
    Dt05 = dt_3_2;
}
double TimeStep::mCalcSoundSpeed(double pressure, double density)
{
    return pow(mGamma * pressure/density, 0.5);
}
double TimeStep::mCourant(double sound_speed, double dx)
{
    return(0.5 * (dx/sound_speed)); 
}
double TimeStep::mVolumePressureTimeLimit(double volume, double pressure, double dx)
{
    return((pow(pressure * volume, 0.5) / dx));
}
double TimeStep::mVolumeTimeLimit(double volumeT, double volumeT_1)
{
    double vol_n_1_2 = (volumeT + volumeT_1) / 2;
    return((volumeT - volumeT_1)/vol_n_1_2);
}
double TimeStep::mRadiationTimeLimit(double RadiationT, double RadiationT_1)
{
    double rad_n_1_2 = (RadiationT + RadiationT_1) / 2; 
    return((RadiationT - RadiationT_1)/rad_n_1_2);
}
double TimeStep::mElectronTemperatureTimeLimit(double electron_temperatureT, double electron_temperatureT_1)
{
    double Te_n_1_2 = (electron_temperatureT + electron_temperatureT_1) / 2; 
    return((electron_temperatureT - electron_temperatureT_1)/Te_n_1_2);
}
double TimeStep::mIonTemperatureTimeLimit(double ion_temperatureT, double ion_temperatureT_1)
{
    double Ti_n_1_2 = (ion_temperatureT + ion_temperatureT_1) / 2; 
    return((ion_temperatureT - ion_temperatureT_1)/Ti_n_1_2);
}

std::string TimeStep::fileNameCreator(std::string path, std::string var, std::string step, std::string extension)
{
    std::string updated_path = path + "/" + var + "_" + step + extension;
//    V path.append("/");
//     path.append(var);
//     path.append("_");
//     path.append(step);
//     path.append(extension);
    return(updated_path);
}

void TimeStep::init(std::string path)
{
    YAML::Node config = YAML::LoadFile(path);
    Dt05 = config["TimeParameters"]["dt"].as<double>();
    Dt1 = config["TimeParameters"]["dt"].as<double>();
    MaxSteps = config["TimeParameters"]["steps"].as<int>();
    Tmax = config["TimeParameters"]["t_max"].as<double>();
    TotalTime = config["TimeParameters"]["t_init"].as<double>();
    dtGlobalMax = config["TimeParameters"]["dt_max"].as<double>();
    dtGlobalMin = config["TimeParameters"]["dt_min"].as<double>();
    mVolumeLimiter = config["TimeParameters"]["VolLimit"].as<float>();
    mPVLimiter = config["TimeParameters"]["PVLimit"].as<float>();
    mRadiationLimiter = config["TimeParameters"]["RadLimit"].as<float>();
    mElectronTemperatureLimiter =config["TimeParameters"]["ElectronLimit"].as<float>();
    mIonTemperatureLimiter = config["TimeParameters"]["IonLimit"].as<float>();
    mCourantCondition = config["TimeParameters"]["cfl"].as<float>();
    mAdaptiveMethod = config["TimeParameters"]["AdapativeMethod"].as<int>();
    InitialDt = Dt05;
}

void TimeStep::writeOut(std::ofstream &variable, std::vector<double> array, int length)
{
    if (!variable)
    {
        std::cout << "file could not be open for writing ! \n";
    }
    for(int i =0; i < length; i++)
    {
        //variable << std::fixed <<  array[i] << "\n";
        variable << std::fixed << std::setprecision(35) <<array[i] << "\n";
    }
    variable.close(); 
}

void TimeStep::dump(IO &ioData)
{
    std::ofstream Variable(fileNameCreator(ioData.dir_paths["TIME"], "TIME",
                        std::to_string(Step), ".txt"));

    Variable << std::fixed << std::setprecision(50) << TotalTime;
    Variable.close();
}
void TimeStep::setTimeParameters(double pvlimit, double vollimit, double radlimit, double telimit,double tilimit, double cfl)
{
    mVolumeLimiter = pvlimit;
    mPVLimiter = vollimit;
    mRadiationLimiter = radlimit;
    mElectronTemperatureLimiter = telimit;
    mIonTemperatureLimiter = tilimit;
    mCourantCondition = cfl;
}