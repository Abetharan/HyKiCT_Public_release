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
#include "HyKiCT/PlasmaParameters.h"

PlasmaParameters::PlasmaParameters(int nx)
{
    mNx = nx;
}

void PlasmaParameters::calculateCoulombLogEE(GridData *gridData, FixedData const &fixedData)
{
    ///Change
    double CoulombLog;
    for(int i = 0; i < mNx; i++)
    {
        double reduced_ne = gridData->NumberDensityE[i] / 1e20;
        double reduced_Te = (gridData->TemperatureE[i] * fixedData.k_to_ev) / 1e3;

        if(gridData->TemperatureE[i] * fixedData.k_to_ev <= 20)
        {
            CoulombLog = 18 - 0.5 * log(reduced_ne) + 1.5 * log(reduced_Te); 
        }
        else
        {
            CoulombLog = 16 - 0.5 * log(reduced_ne) + log(reduced_Te); 
        }
        if(CoulombLog < 2.5)
        {
            CoulombLog = 2.5;
        }
        gridData->CoulombLogEE[i] = CoulombLog;
    }
}
void PlasmaParameters::calculateCoulombLogEI(GridData *gridData, FixedData const &fixedData)
{
    for(int i =0; i < mNx; i++)
    {
        double cond_2 = 10 * pow(gridData->Zbar[i], 2);
        double cond_3 = gridData->TemperatureE[i] * fixedData.k_to_ev;
        if(cond_3 < cond_2)
        {
            gridData->CoulombLogEI[i] = 23 - log(pow(gridData->NumberDensityE[i] * 1e-6, 0.5) * gridData->Zbar[i] * pow(cond_3, -1.5));
        }
        else
        {
            gridData->CoulombLogEI[i] = 24 - log(pow(gridData->NumberDensityE[i] * 1e-6, 0.5) * pow(cond_3, -1.0));

        }
        if(gridData->CoulombLogEI[i] < 2.5)
        {
            gridData->CoulombLogEI[i] = 2.5;
        }
    }
}
void PlasmaParameters::calculateCoulombLogEI(std::vector<double> &CoulombLogEI, 
                                            std::vector<double> const &TemperatureE, 
                                            std::vector<double> const &NumberDensityE,
                                            std::vector<double> const &Zbar, FixedData const &fixedData)
{
   //Ref NRL 2018 
  for(int i =0; i < CoulombLogEI.size(); i++)
    {
        // double cond_1 = gridData->TemperatureI[i] *fixedData.k_to_ev * fixedData.me_mp_ratio * (1/fixedData.Ar[i]);
        double cond_2 = 10 * pow(Zbar[i], 2);
        double cond_3 = TemperatureE[i] * fixedData.k_to_ev;
        if(cond_3 < cond_2)
        {
            CoulombLogEI[i] = 23 - log(pow(NumberDensityE[i] * 1e-6, 0.5) * Zbar[i] * pow(cond_3, -1.5));
        }
        //if((cond_1 < cond_2) && (cond_2 < cond_3))
        else
        {
            CoulombLogEI[i] = 24 - log(pow(NumberDensityE[i] * 1e-6, 0.5) * pow(cond_3, -1.0));
        }
        
        if(CoulombLogEI[i] < 2.5)
        {
            CoulombLogEI[i] = 2.5;
        }
        // else
        // {
        //     gridData->CoulombLogEI[i] = 16 - log(pow(gridData->NumberDensityI[i], 0.5) * pow(gridData->Zbar[i], 2) * pow(gridData->TemperatureI[i] * fixedData.k_to_ev, -1.5) * (1/fixedData.me_mp_ratio));
        // }
    }


}
void PlasmaParameters::setCoulombLog(GridData *gridData, double coloumblog)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->CoulombLogEE[i] = coloumblog;
        gridData->CoulombLogEI[i] = coloumblog;
    }
}
void PlasmaParameters::calculatePlasmaFrequency(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        double plasma_2 = (gridData->NumberDensityE[i] * pow(fixedData.ELECTRON_CHARGE, 2)) / (fixedData.VACUUM_PERMITTIVITY * fixedData.ELECTRON_MASS);
        gridData->PlasmaFrequency[i] =  pow(plasma_2, 0.5);
    }
}
void PlasmaParameters::calculateMinImpactParameter(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        double p_min_1 = (gridData->Zbar[i] * pow(fixedData.ELECTRON_CHARGE, 2)) / (3 * 4 * M_PI * fixedData.VACUUM_PERMITTIVITY * fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureE[i]);
        double p_min_2 = fixedData.HBAR / (2*pow((3 * fixedData.ELECTRON_MASS * fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureE[i]), 0.5));
        gridData->ImpactParameterMinEI[i] = std::max({p_min_1, p_min_2}); 
    }
}
void PlasmaParameters::calculateThermalVelocity(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        double v_th_2 = (fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureE[i]) / fixedData.ELECTRON_MASS;
        double v_th_i = (fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureI[i]) / (fixedData.Ar[i] * fixedData.PROTON_MASS); 
        gridData->ElectronThermalVelocity[i] =  pow(v_th_2, 0.5);
    }
}
void PlasmaParameters::calculateCoulombLogLaser(GridData *gridData, FixedData const &fixedData, double laserWavelength)
{
    //Johnston and Dawson Phys fluids VOl 16 5 1973 
    //Listed in in NRL fomrulary 2019 also however utilises a max(plasmafreq*impact, freq*impact) which reduces calculations to one step instead of two like done here.  
    double frequency = (2 * M_PI * fixedData.SPEED_OF_LIGHT) / laserWavelength;  
    mCoulombLog.resize(mNx, 0);
    for(int i = 0; i < mNx; i++)
    {
        if(gridData->PlasmaFrequency[i] > 0)
        {
            double lambda_1 = gridData->ElectronThermalVelocity[i] / (gridData->PlasmaFrequency[i] * gridData->ImpactParameterMinEI[i]); 
            double lambda_2 = gridData->ElectronThermalVelocity[i] / (frequency * gridData->ImpactParameterMinEI[i]); 
            mCoulombLog[i] = std::log(std::min({lambda_1,lambda_2}));
            if(mCoulombLog[i] < 2.5)
            {
                mCoulombLog[i] = 2.5;
            }
        }
        else
        {
            mCoulombLog[i] = gridData->CoulombLogEI[i];
        }
    }
}
void PlasmaParameters::calculateLeeMoreCoulombLogs(GridData *gridData, FixedData const &fixedData)
{
    //Lee-More 1984 Phys Fluids
    //HYADES for Degen implementation
    //FLASH for Non-Degen Implementation 
    for(int i = 0; i < mNx; i++)
    {
        //Non-Degen
        // double electron_debye =  (gridData->NumberDensityE[i] * fixedData.ELECTRON_CHARGE *fixedData.ELECTRON_CHARGE) / 
        //                         (4 * fixedData.VACUUM_PERMITTIVITY * fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureE[i]);    
        // double ion_Debye = (gridData->NumberDensityI[i] * gridData->Zbar[i] * gridData->Zbar[i] * fixedData.ELECTRON_CHARGE *fixedData.ELECTRON_CHARGE) / 
        //                     (4 * fixedData.VACUUM_PERMITTIVITY * fixedData.BOLTZMANN_CONSTANT * gridData->TemperatureI[i]);
        // double debye_huckel = 1.0 / (electron_debye + ion_Debye);
        // double impact_ratio = debye_huckel / (gridData->ImpactParameterMinEI[i] * gridData->ImpactParameterMinEI[i]);
        
        //electorn Degenerative effects included
        double T_f = 2.8210197e-15 * pow(gridData->NumberDensityE[i], 2.0/3.0);
        double T_b = pow(gridData->TemperatureE[i]*gridData->TemperatureE[i] + T_f*T_f, 0.5);
        double debye = (gridData->NumberDensityE[i] * fixedData.ELECTRON_CHARGE * fixedData.ELECTRON_CHARGE) / (fixedData.VACUUM_PERMITTIVITY * fixedData.BOLTZMANN_CONSTANT);
        double debye_huckel_degen = pow(1.0 / (debye * (1.0 / T_b + gridData->Zbar[i] / gridData->TemperatureI[i])) , 0.5);
        double R_0 = 0.620351 / pow(gridData->NumberDensityI[i], 0.3333333);
        double b_max = std::max({debye_huckel_degen, R_0});
        double e_4 = (fixedData.ELECTRON_CHARGE*fixedData.ELECTRON_CHARGE*fixedData.ELECTRON_CHARGE*fixedData.ELECTRON_CHARGE) /
                        (4 * M_PI * fixedData.VACUUM_PERMITTIVITY * 4 * M_PI * fixedData.VACUUM_PERMITTIVITY);
        double p_min_1 = (gridData->Zbar[i] * pow(fixedData.ELECTRON_CHARGE, 2)) / (3 * 4 * M_PI * fixedData.VACUUM_PERMITTIVITY * fixedData.BOLTZMANN_CONSTANT * T_b);
        double p_min_2 = fixedData.HBAR / (2*pow((3 * fixedData.ELECTRON_MASS * fixedData.BOLTZMANN_CONSTANT * T_b), 0.5));
        double b_min = p_min_1 + p_min_2;
        double impact_ratio = (b_max*b_max)/(b_min* b_min);
        gridData->CoulombLogEI[i] = 0.5 * log(1 + impact_ratio);
        if(gridData->CoulombLogEI[i] < 2.0)
        {
            gridData->CoulombLogEI[i] = 2.0;
        }

    }


}
void PlasmaParameters::calculateCollisionFrequencyEI(GridData *gridData, FixedData const &fixedData, bool laser)
{
    double coulomb;
    for(int i = 0; i < mNx; i++)
    {
        if(laser)
        {
            coulomb = mCoulombLog[i];
        }
        else
        {
            coulomb = gridData->CoulombLogEI[i];
        }
        gridData->CollisionFrequencyEI[i] = 3.633152417604402e-06 * (gridData->NumberDensityE[i] * gridData->Zbar[i] * coulomb) / pow(gridData->TemperatureE[i], 1.5);
    }
}

void PlasmaParameters::calculateCollisionFrequencyEIOverC(GridData *gridData, FixedData const &fixedData, bool laser)
{
    double coulomb;
    for(int i = 0; i < mNx; i++)
    {
        if(laser)
        {
            coulomb = mCoulombLog[i];
        }
        else
        {
            coulomb = gridData->CoulombLogEI[i];
        }
        gridData->CollisionFrequencyEIOverC[i] = 1.211889198895191e-14 * (gridData->NumberDensityE[i] * gridData->Zbar[i] * coulomb) / pow(gridData->TemperatureE[i], 1.5); 
    }
}

