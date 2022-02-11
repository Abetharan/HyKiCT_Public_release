#ifndef GRIDDATA_HPP
#define GRIDDATA_HPP
#include "Data.hpp"
#include "Switches.h"
#include "FixedData.h"
class GridData : public Data
{
    private:    
        void mPrint(std::vector<double> var, std::string name);
        std::vector<std::string> mDumpVars;
        void mInitMulti(std::string path);
        void mInitDivQ(std::string path);
    public:

        //Cell Wall Centered
        std::vector<double>  Velocity, CellWallCoord, HeatFlowE, HeatFlowI, CoupleOperatorSplitHeatFlowE,//REF Note CoupleOperatorSplit etc is ALways SPITZER-HARM HEAT FLOW 
                             SNBCorrection;                                                             // And HeatFlowE ETC is always Q_VFP
                                                                                                        

        //Cell Centered
        //All Source terms, RadFFEmission,RadFFAbsorption, HeatConductionE/I, Exchange, InvBrem have units of WKg^-1 
        //All other units are standard SI  
        
        std::vector<double> Mass, Viscosity, Density, NumberDensityI, NumberDensityE, 
                            TemperatureI, TemperatureE, InternalEnergyI, InternalEnergyE, PressureE, PressureI, 
                            TotalPressure, SpecificHeatE, DpDtE, SpecificHeatI, DpDtI, HeatConductionI, HeatConductionE, CoupleOperatorSplitHeatConductionE,
                            InverseBrem,  Exchange, PlasmaFrequency, TransmittedLaser, PowerAbsorbed, Zbar, 
                            CellCenteredCoord, RadPressure, RadTemperature, RadPdVPower, RadFFAbsorb, RadFFEmission, RadCompPower,RadFicks,
                            RadEnergyDensity, RadFluxLimCoef, RadHeatCapacity, RadPlanckAbsorptionOpacity, RadPlanckEmissionOpacity,
                            RadRossAbsorptionOpacity, RadCouplingCoefficient, RadDiffusionCoefficient, CoulombLogEI, CoulombLogEE, HeatKappaE, HeatKappaI,
                            ElectronThermalVelocity, ImpactParameterMinEI, CollisionFrequencyEI, CollisionFrequencyEIOverC, TotalElectronSource, TotalIonSource,
                            pDvWork, ElectronLeakSource, IonLeakSource,IonThermalVelocity, TmpALlRadEnergy, ALlRadEnergyDensity, IdealSoundSpeed;

        //dx
        std::vector<double> Dx_k_1_2, Dx_k, SNBDqNL;
        std::vector<std::vector<double>> SNBH;
        //Interpolated Te,ne, Z 
        std::vector<double> InterTemperatureE, InterNumberDensityE, InterZbar, InterCoulombLog;
        std::vector<double> LoadInDivQ;        
        std::vector<double> HeatFlowMultiplier;
        std::vector<double> VFPHeatFlow;
        GridData();        
        GridData(int nx, bool init = false);
        void LinearInterpolate();
        void transferData(std::vector<double> lastData, std::vector<double> currentData);
        void setTotalTime();
        void setMemoryCapacity(int nx, bool init);
        void dump(IO &ioData, int step);
        void dumpRadGroup(IO &ioData, int step, std::vector<double> Energy, int group);
        void init(std::string path) override;
        void init(std::string path, bool initRad, bool couple, int coupleMethod);
        void writeOut(std::ofstream &variable, std::vector<double> array, int length) override;
        void terminalPrint() override;
        std::string fileNameCreator(std::string path, std::string var, std::string extension);
        std::string fileNameCreator(std::string path, std::string var, std::string step, std::string extension);
};
#endif
