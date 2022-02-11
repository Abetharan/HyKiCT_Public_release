//   HyKiCT - 1D Lagrangian Radiation-Hydrodyanmics code for testing Coupling with Kinetic codes.
//   Copyright (C) 2020- Abetharan Antony
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNu General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOuT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICuLAR PuRPOSE.  See the
//   GNu General Public License for more details.
//
//   You should have received a copy of the GNu General Public License
//   along with this program.  If not, see <https://www.gnu.org/licenses/>. 
#include "HyKiCT/FluidDynamics.h"
//Constructor for class Hydro. Initilises all parameters as required.
FluidDynamics::FluidDynamics(int nx)
{
    mNx = nx;
}

void FluidDynamics::updateViscosity(GridData *gridDataT_1,  GridData *gridDataT, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        double c1 = 1.0;//0.5;
        double c2 = 1.0;//0.6666666666;
        if(gridDataT_1->Velocity[i + 1] < gridDataT_1->Velocity[i])
        {
            double cs = pow(fixedData.Gamma * gridDataT_1->TotalPressure[i]/gridDataT_1->Density[i], 0.5);
            double modi_gamma = (fixedData.Gamma + 1)/4.0;
            double delta_v = gridDataT_1->Velocity[i + 1] - gridDataT_1->Velocity[i];
            gridDataT->Viscosity[i] = gridDataT_1->Density[i] * abs(delta_v) * 
                                    (c2 * modi_gamma * abs(delta_v) + pow(pow(c2,2)*pow(modi_gamma, 2) * pow(delta_v, 2) + pow(c1*cs,2), 0.5)); //Kuropatenko
                                    
            // gridDataT->Viscosity[i] = fixedData.Cq * gridDataT_1->Density[i] *
            //                          std::pow((gridDataT_1->Velocity[i + 1] - gridDataT_1->Velocity[i]), 2);   /Von-neumann Richtmyer         
            if(gridDataT->Viscosity[i] < mViscosityCap)
            {
                gridDataT->Viscosity[i] = mViscosityCap;
            }
        }
        else
        {
            gridDataT->Viscosity[i] = 0.0;
        }
    }
}

void FluidDynamics::updateMomentum(GridData *gridDataT_1,  GridData *gridDataT, FixedData const &fixedData, TimeStep const &timeData)
{
    double m_k_0_5;
    for(int i = 1; i < mNx; i++)
    {
        m_k_0_5 = 0.5 * (gridDataT_1->Mass[i] + gridDataT_1->Mass[i - 1]);
        gridDataT->Velocity[i] = gridDataT_1->Velocity[i] - (timeData.Dt1/m_k_0_5) * 
                                (gridDataT_1->TotalPressure[i] + gridDataT->Viscosity[i] -
                                 gridDataT->Viscosity[i - 1] - gridDataT_1->TotalPressure[i - 1]);
        if(gridDataT->Velocity[i] <mVelocityCap)
        {
            gridDataT->Velocity[i] = mVelocityCap;
        }
    }

    //******Boundary Conditions*******//
    gridDataT->Velocity[0] = 0;
    
    if(fixedData.RightFluidBoundaryCondition == "r")
    {
        gridDataT->Velocity[mNx] = 0;
    }
    else if(fixedData.RightFluidBoundaryCondition == "f")
    {
        double m_k_0_5 = gridDataT_1->Mass[mNx - 1];
        gridDataT->Velocity[mNx] = gridDataT_1->Velocity[mNx] - (timeData.Dt1/m_k_0_5) * (- gridDataT->Viscosity[mNx - 1] -gridDataT_1->TotalPressure[mNx - 1]);
        if(gridDataT->Velocity[mNx] < 1e-30)
        {
            gridDataT->Velocity[mNx] = 0.0;
        }
    }
}

void FluidDynamics::updateCoords(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData)
{   
    for(int i = 0; i <= mNx; i++)
    {
        gridDataT->CellWallCoord[i] = gridDataT_1->CellWallCoord[i] + gridDataT->Velocity[i] * timeData.Dt05;
        assert(gridDataT->CellWallCoord[i] >= 0);
    }       
}
void FluidDynamics::updateCellCentreCoords(GridData *gridDataT)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->CellCenteredCoord[i] = (gridDataT->CellWallCoord[i + 1] + gridDataT->CellWallCoord[i]) / 2;
    }
}
void FluidDynamics::updateDxs(GridData *gridDataT)
{
    for(int i = 0; i < mNx; i++)
    {
        if(i != mNx - 1)
        {
            gridDataT->Dx_k[i] = gridDataT->CellCenteredCoord[i + 1] - gridDataT->CellCenteredCoord[i]; //delta x_{k} //wall centered
        }
        gridDataT->Dx_k_1_2[i] = gridDataT->CellWallCoord[i + 1] - gridDataT->CellWallCoord[i]; //delta x_{k+1/2} // centered
    }
}

void FluidDynamics::updateDensity(GridData *gridDataT)
{   
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->Density[i] = gridDataT->Mass[i] / gridDataT->Dx_k_1_2[i];
        assert(gridDataT->Density[i] >= 0);      
    }
}

void FluidDynamics::updateNumberDensityI(GridData *gridDataT, FixedData const &fixedData)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->NumberDensityI[i] = gridDataT->Density[i] / (fixedData.Ar[i] * fixedData.PROTON_MASS);
    }
}

void FluidDynamics::updateNumberDensityE(GridData *gridDataT)
{
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->NumberDensityE[i] = gridDataT->NumberDensityI[i] * gridDataT->Zbar[i];
    }
}

void FluidDynamics::updateTemperatureI(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic)
{
    double source = 0;
    double deltaRho;
    double pDvWork;
    for(int i = 0; i < mNx; i++)
    {
        deltaRho = (1/gridDataT->Density[i]) - (1/gridDataT_1->Density[i]);
        pDvWork = (gridDataT_1->TemperatureI[i] * gridDataT_1->DpDtI[i] + gridDataT->Viscosity[i])*deltaRho;
        gridDataT->pDvWork[i] = pDvWork;
        if(!adiabatic)
        {
            source =  gridDataT_1->HeatConductionI[i] - gridDataT_1->Exchange[i];
        }
        gridDataT->TotalIonSource[i] = source;
        gridDataT->TemperatureI[i] =  gridDataT_1->TemperatureI[i] + (1/gridDataT_1->SpecificHeatI[i]) * (timeData.Dt05 * source - pDvWork);
        assert(gridDataT->TemperatureI[i] > 0);
    }
}

void FluidDynamics::updateTemperatureIOperatorSplit(GridData *gridDataT_1, GridData *gridDataT, TimeStep const &timeData)
{
    double deltaRho;
    double pDvWork;
    for(int i = 0; i < mNx; i++)
    {
        deltaRho = (1/gridDataT->Density[i]) - (1/gridDataT_1->Density[i]);
        pDvWork = (gridDataT_1->TemperatureI[i] * gridDataT_1->DpDtI[i] + gridDataT->Viscosity[i])*deltaRho;
        gridDataT->TemperatureI[i] =  gridDataT_1->TemperatureI[i] + (1/gridDataT_1->SpecificHeatI[i]) * (-1* pDvWork);
        assert(gridDataT->TemperatureI[i] > 0);
    }
}
void FluidDynamics::updateTemperatureIOperatorSplit(GridData *gridDataT, TimeStep const &timeData, std::vector<double> *source)
{
    std::vector<double> deref_source = *source;
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->TemperatureI[i] +=  (1/gridDataT->SpecificHeatI[i]) * (timeData.Dt05 * deref_source[i]);
        assert(gridDataT->TemperatureE[i] >= 0);
    }
}    
void FluidDynamics::updateTemperatureE(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic)                   
{
    double source = 0;
    double deltaRho;
    double pDvWork;
    for(int i = 0; i < mNx; i++)
    {
        deltaRho = (1/gridDataT->Density[i]) - (1/gridDataT_1->Density[i]);
        pDvWork = (gridDataT_1->TemperatureE[i] * gridDataT_1->DpDtE[i])*deltaRho;
        // gridDataT->TemperatureE[i] =  gridDataT_1->TemperatureE[i] - (1/gridDataT_1->SpecificHeatE[i]) * pDvWork;
        if(!adiabatic)
        {
            source =  gridDataT_1->HeatConductionE[i] + gridDataT_1->InverseBrem[i] 
                    + gridDataT_1->Exchange[i] + gridDataT_1->RadFFAbsorb[i] - gridDataT_1->RadFFEmission[i];
        } 
        gridDataT->TotalElectronSource[i] = source;
        gridDataT->TemperatureE[i] =  gridDataT_1->TemperatureE[i] + (1/gridDataT_1->SpecificHeatE[i]) * (timeData.Dt05 * source - pDvWork);
        assert(gridDataT->TemperatureE[i] >= 0);
    }
}    
void FluidDynamics::updateTemperatureEOperatorSplit(GridData *gridDataT_1, GridData *gridDataT, TimeStep const &timeData)
{
    double deltaRho;
    double pDvWork;
    for(int i = 0; i < mNx; i++)
    {
        deltaRho = (1/gridDataT->Density[i]) - (1/gridDataT_1->Density[i]);
        pDvWork = (gridDataT_1->TemperatureE[i] * gridDataT_1->DpDtE[i])*deltaRho;
        // gridDataT->TemperatureE[i] =  gridDataT_1->TemperatureE[i] - (1/gridDataT_1->SpecificHeatE[i]) * pDvWork;
        gridDataT->TemperatureE[i] =  gridDataT_1->TemperatureE[i] + (1/gridDataT_1->SpecificHeatE[i]) * (-1 * pDvWork);
        assert(gridDataT->TemperatureE[i] >= 0);
    }
}
void FluidDynamics::updateTemperatureEOperatorSplit(GridData *gridDataT, TimeStep const &timeData, std::vector<double> *source)
{
    std::vector<double> deref_source = *source;
    for(int i = 0; i < mNx; i++)
    {
        gridDataT->TemperatureE[i] += (1/gridDataT->SpecificHeatE[i]) * (timeData.Dt05 * deref_source[i]);
        assert(gridDataT->TemperatureE[i] >= 0);
    }
}    

void FluidDynamics::updateTemperatureCoupleOperatorSplit(GridData *gridDataT_1,  GridData *gridDataT, TimeStep const &timeData, bool adiabatic)
{
    double source = 0;
    for(int i = 0; i < mNx; i++)
    {
        if(!adiabatic)
        {
            source =  gridDataT_1->HeatConductionE[i] - gridDataT_1->CoupleOperatorSplitHeatConductionE[i]; // Note that HeatConduction is the Couple Heat-Conduction and 
        }                                                                                                   // CoupleOperatorSplitHeatConduction is Spitzer-Harm

        gridDataT->TemperatureE[i] =  gridDataT_1->TemperatureE[i] + (1/gridDataT_1->SpecificHeatE[i]) * (timeData.Dt05 * source);
        assert(gridDataT->TemperatureE[i] >= 0);
    }

}
void FluidDynamics::updateTemperatureImplicitHeatConduction(std::vector<double> &Temperature, std::vector<double> &HeatKappa, std::vector<double> const &SpecificHeat, std::vector<double> const &Mass,
                                                           std::vector<double> const &Dx_k, double dt05, std::vector<double> const *multipliers,  bool CoupleMulti)
{
    if(CoupleMulti)
    {
        std::transform(HeatKappa.begin(), HeatKappa.end(), multipliers->begin(), HeatKappa.begin(), std::multiplies<double>());
    }
    double dx_diff, dx_wall_back, dx_wall_forward, source;
    std::vector<double> diffusion_coeficient;
    std::vector<std::tuple<double, double, double>> diff_coef;
    std::tuple<double, double> leftbc, rightbc;
    double val_0,val_1,val_2;
    std::vector<double> source_vec;
    source_vec.resize(mNx);
    diff_coef.resize(mNx);
    diffusion_coeficient.resize(mNx + 1);
    
    double a[mNx];
    double b[mNx];
    double c[mNx]; 
    double d[mNx];
    for(int i = 1; i< mNx; i++)
    {
        diffusion_coeficient[i] = (HeatKappa[i] * dt05) / (SpecificHeat[i] * Dx_k[i - 1] * Mass[i]);
    }
    for(int i = 1; i < mNx - 1; i++)
    {
        val_0 = -1*diffusion_coeficient[i];
        val_1 = diffusion_coeficient[i + 1] + diffusion_coeficient[i] + 1;
        val_2 = -1*diffusion_coeficient[i + 1];
        diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
        source_vec[i] = Temperature[i];
        a[i] = val_0;
        b[i] = val_1;
        c[i] = val_2;
        d[i] = Temperature[i];
    }
    
    int i = 0;
    std::get<0>(leftbc) = 1 + diffusion_coeficient[i + 1];
    std::get<1>(leftbc) = -1*diffusion_coeficient[i + 1];
    source_vec[i] = Temperature[i];
    b[i] = 1 + diffusion_coeficient[i + 1];
    c[i] = -1*diffusion_coeficient[i + 1];
    d[i] = Temperature[i]; 

    i = mNx - 1;
    std::get<0>(rightbc) = -1*diffusion_coeficient[i]; 
    std::get<1>(rightbc) = 1 + diffusion_coeficient[i];
    source_vec[i] = Temperature[i];
    a[i] = -1*diffusion_coeficient[i]; 
    b[i] = 1 + diffusion_coeficient[i];
    d[i] = Temperature[i]; 

    // matrixSolver.fillTriMatrix(diff_coef, source_vec, leftbc, rightbc, false, false);
    // matrixSolver.solveMatrix()
    // matrixSolver.returnSolution(gridDataT->TemperatureE);
    mSolveTridiagonal(a, b, c, d, mNx);
    // //Extract value from PETSC vector and store into std::vector 
    // PetscScalar value_1;
    for(int j = 0; j < mNx; j++)
    {
        Temperature[j] = d[j];//fillVector[j];
        
    }
}
// void FluidDynamics::updateTemperatureImplicitHeatConduction(GridData *gridDataT, TimeStep const &timeData, Switches const &switchData)
// {
//     if(switchData.CoupleMulti)
//     {
//         std::transform(gridDataT->HeatKappaE.begin(), gridDataT->HeatKappaE.end(), gridDataT->HeatFlowMultiplier.begin(), gridDataT->HeatKappaE.begin(), std::multiplies<double>());
//     }
//     double dx_diff, dx_wall_back, dx_wall_forward, source;
//     std::vector<double> diffusion_coeficient;
//     std::vector<std::tuple<double, double, double>> diff_coef;
//     std::tuple<double, double> leftbc, rightbc;
//     double val_0,val_1,val_2;
//     std::vector<double> source_vec;
//     source_vec.resize(mNx);
//     diff_coef.resize(mNx);
//     diffusion_coeficient.resize(mNx);
    
//     double a[mNx];
//     double b[mNx];
//     double c[mNx]; 
//     double d[mNx];
//     for(int i = 1; i< mNx - 1; i++)
//     {
//         diffusion_coeficient[i] = (gridDataT->HeatKappaE[i] * timeData.Dt05) / (gridDataT->SpecificHeatE[i] * gridDataT->Dx_k[i] * gridDataT->Mass[i]);
//     }
//     for(int i = 1; i < mNx - 1; i++)
//     {
//         val_0 = -1*diffusion_coeficient[i];
//         val_1 = diffusion_coeficient[i + 1] + diffusion_coeficient[i] + 1;
//         val_2 = -1*diffusion_coeficient[i + 1];
//         diff_coef[i] = std::make_tuple(val_0, val_1, val_2);
//         source_vec[i] = gridDataT->TemperatureE[i];
//         a[i] = val_0;
//         b[i] = val_1;
//         c[i] = val_2;
//         d[i] = gridDataT->TemperatureE[i];
//     }
    
//     int i = 0;
//     std::get<0>(leftbc) = 1 + diffusion_coeficient[i + 1];
//     std::get<1>(leftbc) = -1*diffusion_coeficient[i + 1];
//     source_vec[i] = gridDataT->TemperatureE[i];
//     b[i] = 1 + diffusion_coeficient[i + 1];
//     c[i] = -1*diffusion_coeficient[i + 1];
//     d[i] = gridDataT->TemperatureE[i]; 

//     i = mNx - 1;
//     std::get<0>(rightbc) = -1*diffusion_coeficient[i]; 
//     std::get<1>(rightbc) = 1 + diffusion_coeficient[i];
//     source_vec[i] = gridDataT->TemperatureE[i];
//     a[i] = -1*diffusion_coeficient[i]; 
//     b[i] = 1 + diffusion_coeficient[i];
    // d[i] = gridDataT->TemperatureE[i]; 

    // matrixSolver.fillTriMatrix(diff_coef, source_vec, leftbc, rightbc, false, false);
    // matrixSolver.solveMatrix()
    // matrixSolver.returnSolution(gridDataT->TemperatureE);
    // mSolveTridiagonal(a, b, c, d, mNx);
    // //Extract value from PETSC vector and store into std::vector 
    // PetscScalar value_1;
    // for(int j = 0; j < mNx; j++)
    // {
        //  gridDataT->TemperatureE[j] = d[j];//fillVector[j];
    // }
// }

void FluidDynamics::mSolveTridiagonal(double* a, double* b, double* c, double* d, int n) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}

void FluidDynamics::fluidLeakSources(GridData *gridData, FixedData const &fixedData, double area, double start_index, double end_index)
{
    double summed_density = 0;
    for(int i = start_index; i < end_index; i++)
    {
        summed_density += (1/gridData->Dx_k_1_2[i]);
    }
    for(int i = start_index; i<end_index; i++)
    {
        gridData->ElectronLeakSource[i] = 0*0.1e3 * area *  gridData->InternalEnergyE[i] * gridData->ElectronThermalVelocity[i] * summed_density; 
        gridData->IonLeakSource[i] = 0*0.1e3 * area * gridData->InternalEnergyI[i] * gridData->IonThermalVelocity[i] * summed_density;
    }
}