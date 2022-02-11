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
#include "HyKiCT/IntSource.h"

IntSource::IntSource(int nx)
{
    mNx = nx;
}
double IntSource::mCalculateKappaE(double TemperatureE, double CoulombLogEI, double Zbar)
{

    // double ke = 1.843076667547614E-10 * std::pow(TemperatureE, 2.5) * std::pow(CoulombLogEI, -1)
    //                                 *  std::pow(Zbar, -1); // Checked 27/01/2020 only true for Z = 1
    // double ke =  13.6*((Zbar+0.24)/(Zbar+4.2))*5.759614586E-11* pow(TemperatureE, 2.5) * pow(CoulombLogEI, -1) *  pow(Zbar, -1);
    double ke =  13.6*((Zbar+0.24)/(Zbar+4.2))*5.759614586E-11*  TemperatureE * TemperatureE * sqrt(TemperatureE) * (1.0/CoulombLogEI) * (1.0 / Zbar);
    return ke;
}
void IntSource::calculateLeeMoreKappaECorr(GridData *gridData, FixedData const &fixedData)
{
    //Lee-more 1984 Phys Fluids 
    //Based on Hyades implementation
    //We skip corrections the collision time i.e the ln(1+exp(mukt)F_1/2 term.
    //This adds high lvl of complexity for mininimal correction value .. dont suspect this will 
    //more detail. The coloumb log corrections prob the most important. 
    for(int i = 0; i < mNx; i++)
    {
        double T_f = 2.8210197e-15 * pow(gridData->NumberDensityE[i], 2.0/3.0);
        double x = pow(T_f / gridData->TemperatureE[i], 0.5);
        double degen_param = x*x*x*((0.7531 + 0.1679 * x + 0.3108 * x * x) / ( 1+ 0.2676 * x + 0.2280 * x * x + 0.3099 * x * x* x));
        double A_beta = (13.5 + degen_param * (0.976 + 0.436 * degen_param)) / (1 + degen_param * (0.51 + 0.126 * degen_param));
        gridData->HeatKappaE[i] *= A_beta/13.6;  //Dividing by 13.6 as A-beta incompasse the standard Spitzer/E-H coeficient 
    }
}
void IntSource::calculateImplicitKappaE(GridData *gridData, FixedData const &fixedData)
{
    //For readabilities sake the un split version is 
    //fixedData.heatflxlmE * InvSqrt(fixedData.ELECTRON_MASS) * pow(fixedData.BOLTZMANN_CONSTANT, 1.5)  
    //                                * gridData->NumberDensityE[i] * pow(gridData->TemperatureE[i], 1.5);
    for(int i = 1; i < mNx; i++)
    {
;
        double centered_ke = 2 * (gridData->HeatKappaE[i] * gridData->HeatKappaE[i - 1]) / (gridData->HeatKappaE[i] + gridData->HeatKappaE[i - 1]);
        // centered_ke /= (1.0 + 1.0/(gridData->CoulombLogEI[i] *6));
        double Te_flux =  ((gridData->TemperatureE[i] - gridData->TemperatureE[i - 1]) /
                                 (gridData->CellCenteredCoord[i] - gridData->CellCenteredCoord[i - 1]));
        double e_flux  = centered_ke * abs(Te_flux);

        double Te_boundary = 0.5 * (gridData->TemperatureE[i] + gridData->TemperatureE[i - 1]);
        double ne_boundary = 0.5 * (gridData->NumberDensityE[i] + gridData->NumberDensityE[i - 1]);
        double e_flux_lim =  fixedData.heatflxlmE * InvSqrt(fixedData.ELECTRON_MASS)  * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
                                    * ne_boundary * Te_boundary*sqrt(Te_boundary);
        // double e_flux_lim =  fixedData.heatflxlmE * InvSqrt(fixedData.ELECTRON_MASS) * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
        //                             * gridData->NumberDensityE[i] * gridData->TemperatureE[i]*sqrt(gridData->TemperatureE[i]);
        
        if(fixedData.heatflxlmE < 1.0)
        {
            double ke_eff = centered_ke / ( 1 + e_flux / e_flux_lim); 
            gridData->HeatKappaE[i] = ke_eff;
        }
        else
        {
            gridData->HeatKappaE[i] = centered_ke;
        }
    }
}
void IntSource::calculateImplicitKappaI(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 1; i < mNx; i++)
    {
        double centered_ki = 2 * (gridData->HeatKappaI[i] * gridData->HeatKappaI[i - 1]) / (gridData->HeatKappaI[i] + gridData->HeatKappaI[i - 1]);
        double Ti_flux =  ((gridData->TemperatureI[i] - gridData->TemperatureI[i - 1])
                                                    / (gridData->CellCenteredCoord[i] - gridData->CellCenteredCoord[i - 1]));
        double i_flux = centered_ki * Ti_flux;
        if((fixedData.heatflxlmI < 1.0))
        {
            
            // double i_flux_lim =  fixedData.heatflxlmI * InvSqrt(fixedData.Ar[i]*fixedData.PROTON_MASS) * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
            //                         * gridData->NumberDensityI[i] * gridData->TemperatureI[i]*sqrt(gridData->TemperatureI[i]);;
            
            double Ti_boundary = 0.5 * (gridData->TemperatureI[i] + gridData->TemperatureI[i - 1]);
            double ni_boundary = 0.5 * (gridData->NumberDensityI[i] + gridData->NumberDensityI[i - 1]);
            double Ar_boundary = 0.5 * (fixedData.Ar[i] + fixedData.Ar[i - 1]);
            double i_flux_lim =  fixedData.heatflxlmI * InvSqrt(Ar_boundary*fixedData.PROTON_MASS) * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
                                    * ni_boundary * Ti_boundary*sqrt(Ti_boundary);
            double ki_eff = centered_ki / ( 1 + i_flux / i_flux_lim); 
            gridData->HeatKappaI[i] =  ki_eff;
        }
        else
        {
            gridData->HeatKappaI[i] = centered_ki;
        }
    }
}


double IntSource::mCalculateKappaI(double TemperatureI, double CoulombLogEI, double Zbar, int Ar)
{
    //Readable form
    // double ki = 5.2420798880426694e-12* std::pow(TemperatureI, 2.5) * (1.0 / CoulombLogEI) 
    //                         * InvSqrt(Ar) * std::pow(Zbar, -4);
    double ki = 5.2420798880426694e-12* TemperatureI * TemperatureI *sqrt(TemperatureI) * (1.0 / CoulombLogEI) 
                            * InvSqrt(Ar) * (1.0/(Zbar*Zbar*Zbar*Zbar));
    return ki;
}
void IntSource::calculateKappa(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridData->HeatKappaE[i] = IntSource::mCalculateKappaE(gridData->TemperatureE.at(i), gridData->CoulombLogEI.at(i), gridData->Zbar.at(i));
        gridData->HeatKappaI[i] = IntSource::mCalculateKappaI(gridData->TemperatureI.at(i), gridData->CoulombLogEI.at(i), gridData->Zbar.at(i), fixedData.Ar.at(i));
    }
}


void IntSource::diagonsticHeatFlowE(GridData *gridData, FixedData const &fixedData, Switches const &switchData)
{
    for(int i = 1; i< mNx; i++)
    {
        double centered_ke = gridData->HeatKappaE[i];
        gridData->HeatFlowE[i] = centered_ke * ((gridData->TemperatureE[i] - gridData->TemperatureE[i - 1]) /
                                 (gridData->CellCenteredCoord[i] - gridData->CellCenteredCoord[i - 1]));
        if(switchData.CoupleMulti)
        {
            gridData->HeatFlowE[i] *= gridData->HeatFlowMultiplier[i];
        }
    }

}
void IntSource::heatFlowE(GridData *gridData, FixedData const &fixedData, Switches const &switchData)
{
    //Note the imposition of no heat flow from boundaries are imposed in the initial creation of gridData
    double e_flux_lim;
    for(int i = 1; i< mNx; i++)
    {
        double centered_ke = 2 * (gridData->HeatKappaE[i] * gridData->HeatKappaE[i - 1]) / (gridData->HeatKappaE[i] + gridData->HeatKappaE[i - 1]);
        double Te_flux =  ((gridData->TemperatureE[i] - gridData->TemperatureE[i - 1]) /
                                 (gridData->CellCenteredCoord[i] - gridData->CellCenteredCoord[i - 1]));
        double e_flux  = centered_ke * Te_flux;

        //Readable Form:
        // e_flux_lim =  fixedData.heatflxlmE * InvSqrt(fixedData.ELECTRON_MASS) * std::pow(fixedData.BOLTZMANN_CONSTANT, 1.5)  
        //                             * gridData->NumberDensityE[i] * std::pow(gridData->TemperatureE[i], 1.5);
        double Te_boundary = 0.5 * (gridData->TemperatureE[i] + gridData->TemperatureE[i - 1]);
        double ne_boundary = 0.5 * (gridData->NumberDensityE[i] + gridData->NumberDensityE[i - 1]);
        e_flux_lim =  fixedData.heatflxlmE * InvSqrt(fixedData.ELECTRON_MASS)  * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
                                    * ne_boundary * Te_boundary*sqrt(Te_boundary);
        if((switchData.SNBHeatFlow)||(switchData.Couple))
        {
            gridData->HeatFlowE[i] = e_flux;
        }
        else
        {   
            if((fixedData.heatflxlmE < 1.0) && (fixedData.heatflxlmE > 0.0)) 
            {
                double ke_eff = centered_ke / ( 1 + abs(e_flux) / e_flux_lim); // Note that the explicit version does not abs(Te_flux) like the implicit. 
                gridData->HeatFlowE[i] = ke_eff * Te_flux; 
            }
            else
            {
                gridData->HeatFlowE[i] = e_flux;
            }
        }
    }
}

void IntSource::heatFlowI(GridData *gridData, FixedData const &fixedData)
{
    for(int i = 1; i< mNx; i++)
    {
        // double centered_ki = 0.5 * (gridData->HeatKappaI[i] + gridData->HeatKappaI[i - 1]);
        double centered_ki = 2 * (gridData->HeatKappaI[i] * gridData->HeatKappaI[i - 1]) / (gridData->HeatKappaI[i] + gridData->HeatKappaI[i - 1]);
        double Ti_flux =  ((gridData->TemperatureI[i] - gridData->TemperatureI[i - 1])
                                                    / (gridData->CellCenteredCoord[i] - gridData->CellCenteredCoord[i - 1]));
        double i_flux = centered_ki * Ti_flux;
        if((fixedData.heatflxlmI < 1.0))
        {
            
            // double i_flux_lim =  fixedData.heatflxlmI * InvSqrt(fixedData.Ar[i]*fixedData.PROTON_MASS) * std::pow(fixedData.BOLTZMANN_CONSTANT, 1.5)  
            //                             * gridData->NumberDensityI[i] * std::pow(gridData->TemperatureI[i], 1.5); //Not consistent. Use cell-Walls as below
            double Ti_boundary = 0.5 * (gridData->TemperatureI[i] + gridData->TemperatureI[i - 1]);
            double ni_boundary = 0.5 * (gridData->NumberDensityI[i] + gridData->NumberDensityI[i - 1]);
            double Ar_boundary = 0.5 * (fixedData.Ar[i] + fixedData.Ar[i - 1]);
            double i_flux_lim =  fixedData.heatflxlmI * InvSqrt(Ar_boundary*fixedData.PROTON_MASS) * fixedData.BOLTZMANN_CONSTANT * sqrt(fixedData.BOLTZMANN_CONSTANT)  
                                    * ni_boundary * Ti_boundary*sqrt(Ti_boundary);
            
            double ki_eff = centered_ki / ( 1.0 + i_flux / i_flux_lim); 
            gridData->HeatFlowI[i] =  ki_eff * Ti_flux;//std::min({i_flux, i_flux_lim});
        }
        else
        {
            gridData->HeatFlowI[i] = i_flux;
        }
    }
}
void IntSource::subtractHeatFlowE(GridData *gridData)
{
    for(int i = 0; i<=mNx; i++)
    {
        gridData->HeatFlowE[i] = gridData->HeatFlowE[i]  - gridData->VFPHeatFlow[i];
    }

}
void IntSource::multiplierHeatFlowE(GridData *gridData, FixedData const&fixedData)
{
  //heatFlowE(gridData, fixedData);
    for(int i = 0; i<=mNx; i++)
    {
        gridData->HeatFlowE[i]*=gridData->HeatFlowMultiplier[i];
    }
}
void IntSource::thermalConducE(GridData *gridData)
{
    for(int i = 0; i< mNx; i ++)
    {
        gridData->HeatConductionE[i] = (gridData->HeatFlowE[i + 1] - gridData->HeatFlowE[i]) / gridData->Mass[i];
    }
}

void IntSource::operatorSplitThermalConducE(GridData *gridData)
{
    for(int i = 0; i< mNx; i ++)
    {
        gridData->CoupleOperatorSplitHeatConductionE[i] = (gridData->CoupleOperatorSplitHeatFlowE[i + 1] - gridData->CoupleOperatorSplitHeatFlowE[i]) / gridData->Mass[i];
    }
}
void IntSource::thermalConducI(GridData *gridData)
{
    for(int i = 0; i< mNx; i ++)
    {
        gridData->HeatConductionI[i] = (gridData->HeatFlowI[i + 1] - gridData->HeatFlowI[i]) / gridData->Mass[i];
    }
}

// void IntSource::exchangeCoefficient(GridData *gridData, FixedData const &fixedData, TimeStep const &timeData)
// {
    // for(int i = 0; i < mNx; i++)
    // {
        // standard_Exchange = 4.9044757184917434e-05 * std::pow(gridData->Zbar[i], 2) * std::pow(fixedData.Ar[i], -1) * gridData->CoulombLogEI[i] * gridData->NumberDensityE[i] * 
        //                     std::pow(gridData->TemperatureE[i], -1.5);         
    // }

// }
void IntSource::exchange(GridData *gridData, FixedData const &fixedData, TimeStep const&timeData)
{
    // Solving 
    // Ref: NRL formulary 2019, HYADES 
    // 3(me/mp)n_ek_b / tau_ei (Ti - Te)  
    // Notice the 1/density term included to make the units WKg^_1 as required. 
    // We cap the exchange term to be change max energy change allowed by specific heat in a single time step.
    // This is needed at very small temperatures where the exchange term can be untractably small for a explicit method.
    // Ref: GORGON. 
    double standard_Exchange;
    double dute_exchange;
    double duti_exchange;
    float qf = 0.5;
    for(int i = 0; i < mNx; i++)
    {
        //Readable Form
        // standard_Exchange = 4.9044757184917434e-05 * std::pow(gridData->Zbar[i], 2) * (1.0 / fixedData.Ar[i]) * gridData->CoulombLogEI[i] * gridData->NumberDensityE[i] * 
        //                     std::pow(gridData->TemperatureE[i], -1.5);         
        // Equivalent formulation to the one being used. Here the option to use collision frequency directly is included 
        // This will allow for the complete Lee-More implementation if that is ever needed i.e. adding Fermi integrals to the collision
        // time to account for degenerate effects. At this point only included through Lee-More Coulomb Logs. 
        // standard_Exchange = 3 * (1.0 / gridData->Density[i]) * fixedData.BOLTZMANN_CONSTANT* fixedData.me_mp_ratio * (1.0 / fixedData.Ar[i]) * gridData->NumberDensityE[i]  * gridData->CollisionFrequencyEI[i];
        standard_Exchange = 4.9044757184917434e-05 * gridData->Zbar[i] * gridData->Zbar[i] * (1.0 /(fixedData.Ar[i] * fixedData.Ar[i])) * gridData->CoulombLogEI[i] * gridData->NumberDensityE[i] * 
                            InvSqrt(gridData->TemperatureE[i]) * (1.0/gridData->TemperatureE[i]);       
        dute_exchange = (qf * gridData->SpecificHeatE[i])  /  (2 * timeData.Dt05);
        duti_exchange = (qf * gridData->SpecificHeatI[i])  / (2 * timeData.Dt05);
        gridData->Exchange[i] = (gridData->TemperatureI[i] - gridData->TemperatureE[i]) * std::min({standard_Exchange, dute_exchange, duti_exchange});
        // gridData->Exchange[i] = (gridData->TemperatureI[i] - gridData->TemperatureE[i]) * standard_Exchange;
    }
}
void IntSource::setExchangeCoefficient(GridData *gridData, double ex)
{
    gridData->Exchange[0] = (gridData->TemperatureI[0] -  gridData->TemperatureE[0]) * ex;    
}
