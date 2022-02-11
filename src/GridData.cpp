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
#include "HyKiCT/GridData.h"

GridData::GridData(int nx, bool init)
{
    setMemoryCapacity(nx, init);
}

void GridData::setMemoryCapacity(int nx, bool init)
{
        if(!init)
        {
            Velocity.resize(nx + 1);
            CellWallCoord.resize(nx + 1);
            Mass.resize(nx);
            TemperatureI.resize(nx);
            TemperatureE.resize(nx);
            Density.resize(nx);
            LoadInDivQ.resize(nx);        
            VFPHeatFlow.resize(nx - 1);        
            HeatFlowMultiplier.resize(nx - 1);
        }
        //Cell Wall Centered
        HeatFlowE.resize(nx + 1);
        SNBCorrection.resize(nx + 1);
        CoupleOperatorSplitHeatFlowE.resize(nx + 1);
        HeatFlowI.resize(nx + 1);
        RadDiffusionCoefficient.resize(nx + 1);
        TransmittedLaser.resize(nx + 1);
        SNBDqNL.resize(nx - 1);
        Dx_k.resize(nx - 1);

        //Cell Centered 
        HeatKappaE.resize(nx);
        HeatKappaI.resize(nx);
        CellCenteredCoord.resize(nx);
        Dx_k_1_2.resize(nx);
        Viscosity.resize(nx);
        NumberDensityI.resize(nx);
        NumberDensityE.resize(nx);
        InternalEnergyI.resize(nx);
        InternalEnergyE.resize(nx);
        PressureE.resize(nx);
        PressureI.resize(nx);
        TotalPressure.resize(nx);
        SpecificHeatE.resize(nx);
        SpecificHeatI.resize(nx);
        DpDtE.resize(nx);
        DpDtI.resize(nx);
        HeatConductionI.resize(nx);
        HeatConductionE.resize(nx);
        CoupleOperatorSplitHeatConductionE.resize(nx);
        InverseBrem.resize(nx);
        Exchange.resize(nx);
        PlasmaFrequency.resize(nx);
        ImpactParameterMinEI.resize(nx);
        ElectronThermalVelocity.resize(nx);
        PowerAbsorbed.resize(nx);
        Zbar.resize(nx);
        RadFFEmission.resize(nx);
        RadFFAbsorb.resize(nx);
        RadFicks.resize(nx);
        RadEnergyDensity.resize(nx, 1e-15); 
        RadTemperature.resize(nx, 1E-15);
        CoulombLogEE.resize(nx);
        CoulombLogEI.resize(nx);
        CollisionFrequencyEI.resize(nx);
        CollisionFrequencyEIOverC.resize(nx);
        TotalIonSource.resize(nx);
        TotalElectronSource.resize(nx);
        pDvWork.resize(nx);
        IonLeakSource.resize(nx);
        ElectronLeakSource.resize(nx);
        IonThermalVelocity.resize(nx);
}
std::string GridData::fileNameCreator(std::string path, std::string var, std::string step, std::string extension)
{
    std::string updated_path = path + "/" + var + "_" + step + extension;
//    V path.append("/");
//     path.append(var);
//     path.append("_");
//     path.append(step);
//     path.append(extension);
    return(updated_path);
}

std::string GridData::fileNameCreator(std::string path, std::string var, std::string extension)
{
    std::string updated_path = path + "/" + var + extension;
    return(updated_path);
}

void GridData::init(std::string path, bool initRad, bool couple, int coupleMethod)
{
    GridData::init(path);
    if(initRad)
    {
        std::string file_path = fileNameCreator(path, "rad_energy_density", ".txt");
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            int i = 0;
            while (getline(initFile, strInput))
            {
                TmpALlRadEnergy.push_back(stod(strInput));
                i++;
            }
        }
    }
    if(couple)
    {
        std::string file_path = fileNameCreator(path, "qe", ".txt");
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                switch(coupleMethod){

                    case 1:
                        HeatFlowMultiplier.push_back(stod(strInput));
                        break;
                    case 2:
                        LoadInDivQ.push_back(stod(strInput));
                        break;
                    case 3:
                        VFPHeatFlow.push_back(stod(strInput));
                        break;
                }
            }
        }
    }
}
void GridData::init(std::string path)
{
    std::vector<std::string> vars;
    vars = {"coord", "velocity", "density", "electron_temperature",
                "ion_temperature", "mass"};
    int no_var = vars.size();
    for (int i = 0; i < no_var; i++)
    {
        std::string line;
        std::string file_path = fileNameCreator(path, vars[i], ".txt");
        std::cout << file_path << "\n";
        std::ifstream initFile(file_path);
        std::string strInput;

        if (initFile.is_open())
        {
            while (getline(initFile, strInput))
            {
                if (i == 0)
                {
                    CellWallCoord.push_back(stod(strInput));
                }
                if (i == 1)
                {
                    Velocity.push_back(stod(strInput));
                }
                if (i == 2)
                {
                    Density.push_back(stod(strInput));
                }
                if (i == 3)
                {
                    TemperatureE.push_back(stod(strInput));
                }
                if (i == 4)
                {
                    TemperatureI.push_back(stod(strInput));
                }
                if (i == 5)
                {
                    Mass.push_back(stod(strInput));
                }
               }
        }
    }
}
void GridData::writeOut(std::ofstream &variable, std::vector<double> array, int length)
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

void GridData::dump(IO &ioData, int step)
{
    int total_dump_count = ioData.dumpVars.size();
    std::string dump_variable;
    for(int i = 0; i < total_dump_count; i++)
    {
        if(ioData.dumpBools[ioData.dumpVars[i]])
        {
            dump_variable = ioData.dumpVars[i];
            std::string txt_path = fileNameCreator(ioData.dir_paths[dump_variable], dump_variable, 
                                std::to_string(step), ".txt");
            std::ofstream Variable(txt_path);
            
            if (dump_variable == "VISCOSITY")
            {
                writeOut(Variable, Viscosity, Viscosity.size()); 
            }
            if (dump_variable == "VELOCITY")
            {
                writeOut(Variable, Velocity, Velocity.size()); 
            }
            if (dump_variable == "CELL_WALL_X")
            {
                writeOut(Variable, CellWallCoord, CellWallCoord.size()); 
            }
            if (dump_variable == "CELL_CENTRE_X")
            {
                writeOut(Variable, CellCenteredCoord, CellCenteredCoord.size()); 
            }
            if (dump_variable == "DENSITY")
            {
                writeOut(Variable, Density, Density.size()); 
            }
            if (dump_variable == "ION_NUMBER_DENSITY")
            {
                writeOut(Variable, NumberDensityI, NumberDensityI.size()); 
            }
            if (dump_variable == "ELECTRON_NUMBER_DENSITY")
            {
                writeOut(Variable, NumberDensityE, NumberDensityE.size()); 
            }
            if (dump_variable == "ION_TEMPERATURE")
            {
                writeOut(Variable, TemperatureI, TemperatureI.size()); 
            }
            if (dump_variable == "ELECTRON_TEMPERATURE")
            {
                writeOut(Variable, TemperatureE, TemperatureE.size()); 
            }
            if (dump_variable == "ION_INTERNAL_ENERGY")
            {
                writeOut(Variable, InternalEnergyI, InternalEnergyI.size()); 
            }
            if (dump_variable == "ELECTRON_INTERNAL_ENERGY")
            {
                writeOut(Variable, InternalEnergyE, InternalEnergyE.size()); 
            }
            if (dump_variable == "ELECTRON_PRESSURE")
            {
                writeOut(Variable, PressureE, PressureE.size()); 
            }
            if (dump_variable == "ION_PRESSURE")
            {
                writeOut(Variable, PressureI, PressureI.size()); 
            }
            if (dump_variable == "TOTAL_PRESSURE")
            {
                writeOut(Variable, TotalPressure, TotalPressure.size()); 
            }
            if (dump_variable == "ELECTRON_SPECIFIC_HEAT")
            {
                writeOut(Variable, SpecificHeatE, SpecificHeatE.size()); 
            }
            if (dump_variable == "ELECTRON_DP_DT")
            {
                writeOut(Variable, DpDtE, DpDtE.size()); 
            }
            if (dump_variable == "ION_SPECIFIC_HEAT")
            {
                writeOut(Variable, SpecificHeatI, SpecificHeatI.size()); 
            }
            if (dump_variable == "ION_DP_DT")
            {
                writeOut(Variable, DpDtI, DpDtI.size()); 
            }
            if (dump_variable == "ELECTRON_HEAT_FLOW_X")
            {
                writeOut(Variable, HeatFlowE, HeatFlowE.size()); 
            }
            if (dump_variable == "ION_HEAT_FLOW_X")
            {
                writeOut(Variable, HeatFlowI, HeatFlowI.size()); 
            }
            if (dump_variable == "ELECTRON_HEAT_CONDUCTION")
            {
                writeOut(Variable, HeatConductionE, HeatConductionE.size()); 
            }
            if (dump_variable == "ION_HEAT_CONDUCTION")
            {
                writeOut(Variable, HeatConductionI, HeatConductionI.size()); 
            }
            if (dump_variable == "INVERSE_BREM")
            {
                writeOut(Variable, InverseBrem, InverseBrem.size()); 
            }
            if (dump_variable == "EXCHANGE")
            {
                writeOut(Variable, Exchange, Exchange.size()); 
            }
            if (dump_variable == "PLASMA_FREQUENCY")
            {            
                writeOut(Variable, PlasmaFrequency, PlasmaFrequency.size()); 
            }
            if (dump_variable == "COULOMB_LOG")
            {            
                writeOut(Variable, CoulombLogEI, CoulombLogEI.size()); 
            }
            if (dump_variable == "LASER_STUFF")
            {
                writeOut(Variable, TransmittedLaser, TransmittedLaser.size()); 
            }
            if (dump_variable == "ZBAR")
            {
                writeOut(Variable, Zbar, Zbar.size()); 
            }
            if (dump_variable == "MASS")
            {
                writeOut(Variable, Mass, Mass.size()); 
            }
            if (dump_variable  == "RAD_PRESSURE")
            {
                writeOut(Variable, RadPressure, RadPressure.size()); 
            }
            if (dump_variable  == "RAD_TEMPERATURE")
            {
                writeOut(Variable, RadTemperature, RadTemperature.size()); 
            }
            if (dump_variable  == "RAD_FREE_FREE_EMISSION")
            {
                writeOut(Variable, RadFFEmission, RadFFEmission.size()); 
            }
            if (dump_variable  == "RAD_FREE_FREE_ABSORB")
            {
                writeOut(Variable, RadFFAbsorb, RadFFAbsorb.size()); 
            }
            if (dump_variable  == "RAD_COMP")
            {
                writeOut(Variable, RadCompPower, RadCompPower.size()); 
            }
            if (dump_variable  == "RAD_ENERGY_DENSITY")
            {
                writeOut(Variable, RadEnergyDensity, RadEnergyDensity.size()); 
            }
            if (dump_variable  == "RAD_DIFFUSION")
            {
                writeOut(Variable, RadFicks, RadFicks.size()); 
            }
            if (dump_variable  == "TOTAL_ION_SOURCE")
            {
                writeOut(Variable, TotalIonSource, TotalIonSource.size()); 
            }
            if (dump_variable  == "TOTAL_ELECTRON_SOURCE")
            {
                writeOut(Variable, TotalElectronSource, TotalElectronSource.size()); 
            }
            if (dump_variable  == "PDVWORK")
            {
                writeOut(Variable, pDvWork, pDvWork.size()); 
            }
            if (dump_variable  == "PLANCKOPACITIES")
            {
                writeOut(Variable, RadPlanckAbsorptionOpacity, RadPlanckAbsorptionOpacity.size()); 
            }
            if (dump_variable  == "ROSSLANDOPACITIES")
            {
                writeOut(Variable, RadRossAbsorptionOpacity, RadRossAbsorptionOpacity.size()); 
            }
            if(dump_variable == "RAD_MULTI_GROUP")
            {
                writeOut(Variable, ALlRadEnergyDensity, ALlRadEnergyDensity.size()); 
            }
        }
    }
}
void GridData::dumpRadGroup(IO &ioData, int step, std::vector<double> Energy, int group)
{
    std::string dump_variable = "RAD_ENERGY_DENSITY";
    std::string group_name = "RAD_ENERGY_DENSITY_" + std::to_string(group);
    std::string txt_path = fileNameCreator(ioData.dir_paths[dump_variable], group_name, 
                        std::to_string(step), ".txt");
    std::ofstream Variable(txt_path);
    writeOut(Variable, Energy, Energy.size()); 
}
void GridData::mPrint(std::vector<double> var, std::string name)
{
    //std::cout << std::fixed;
    std::cout << "*--------------------------------------------------------------------------------------*" << "\n";
    std::cout << name <<  " Values : " << " \n" ;
    for(unsigned long int i = 0; i < var.size(); i++ )
    {
        std::cout << var[i] << " ";
    }
    std::cout << "\n*--------------------------------------------------------------------------------------*" << "\n";

}

void GridData::terminalPrint()
{
    std::cout << "*--------------------------------------------------------------------------------------*" << "\n";
    std::cout << "**************************PRINTING ALL HYDRO VARIABLES**********************************" << "\n";
    std::cout << "*--------------------------------------------------------------------------------------*" << std::endl;
       
    mPrint(CellWallCoord, "Cell Wall Coord (m)");
    mPrint(Viscosity, "Viscosity");
    mPrint(Velocity, "Velocity (m/s)");
    mPrint(Mass, "Mass (Kg/m2)");
    mPrint(Density, "Density (Kg/m^-3)");
    mPrint(NumberDensityE, "Number Density E (m^-3)");
    mPrint(NumberDensityI, "Number Density I (m^-3)");
    mPrint(TemperatureE, "Temperature E (K)");
    mPrint(TemperatureI, "Temperature I (K)");
    mPrint(InternalEnergyE, "Internal Energy E (J/Kg)");
    mPrint(InternalEnergyI, "Internal Energy I (J/Kg)");
    mPrint(PressureE, "Pressure E (Pa)");
    mPrint(PressureI, "Pressure I (Pa)");
    mPrint(TotalPressure, "Total Pressure (Pa)");
    mPrint(HeatConductionE, "Heat Conduction E (W/kg)");
    mPrint(HeatConductionI, "Heat Conduction I (W/kg)");
    mPrint(InverseBrem, "Inverse Bremsstrahlun (W/kg)");
    mPrint(Exchange, "Collision Heating(Exchange) (W/Kg)"); 
    mPrint(Zbar, "Ionisation"); 
    mPrint(RadTemperature, "Radiation Temperature");
    mPrint(RadEnergyDensity, "Radiation Energy Density");
    mPrint(RadFFAbsorb, "Radiatio Free-Free-Absorbed");
    mPrint(RadFFEmission, "Radiatio Free-Free-Emission");
    mPrint(RadPdVPower, "Radiatio PdV work");
}

void GridData::LinearInterpolate()
{
    int nx = TemperatureE.size();
    int total_points = 2*nx - 1;
    int j = 1;
    
    if (InterZbar.size() != total_points)
    {
        InterZbar.resize(total_points);
        InterTemperatureE.resize(total_points);
        InterNumberDensityE.resize(total_points);
        InterCoulombLog.resize(total_points);
    }

    for(int i = 0; i < total_points; i++)
    {
        if(i % 2 == 0)
        {
            int index = (i / 2);
            InterTemperatureE[i] = TemperatureE[index] ; 
            InterNumberDensityE[i] = NumberDensityE[index]; 
            InterZbar[i] = Zbar[index];
        }
        else
        {
            double weight = (CellWallCoord[j] - CellCenteredCoord[j - 1]) / 
                        (CellCenteredCoord[j] - CellCenteredCoord[j - 1]);
            InterTemperatureE[i] = TemperatureE[j - 1] * (1 - weight) + TemperatureE[j] * weight; 
            InterNumberDensityE[i] = NumberDensityE[j - 1] * (1 - weight) + NumberDensityE[j] * weight; 
            InterZbar[i] = Zbar[j - 1] * (1 - weight) + Zbar[j] * weight;
            j++;
        }
    }
}