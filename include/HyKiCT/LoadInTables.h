#ifndef LOADINTABLES
#define LOADINTABLES
#include <utility>
#include <functional> 
#include <vector>
#include <array>
#include "GridData.h"
#include "FixedData.h"

class LoadInTables
{
    struct searchIndicies
    {
        std::pair<int, int> te_index;
        std::pair<int, int> rho_index;
        std::pair<double, double> te_values;
        std::pair<double, double> rho_values;
        std::array<int, 4> combined_index;
    };
    private:

            int mNx, mNcols;
            int mStartIndex;
            searchIndicies mIndexSearch(std::vector<double> &temperature, std::vector<double> &density, double search_temperature, double search_density);
            void mTextFileLoader(std::string path, std::vector<double> &mVector, char delimiter = ',');      
            void mFEOSFormatFileLoader(std::string path, std::vector<double> &mVector, bool loadAll = false, char delimiter = ',');
    public:

            LoadInTables(int StartIndex, int nx);
            // LoadInTables(const LoadIntables &Tables) 
            void converTemperatureToeV(double ev_to_conversion);
            double getGradientQuantity(double lowerTemp, double upperTemp, double lowerDensity, double upperDensity, double q11, double q12, double q21, 
                                    double q22, double density);
            void feosEoSLoadIn(std::string path);
            void OpacityLoadIn(std::string path);
            void PhotonGridLoadIn(std::string path);
            void feosIonizationLoadIn(std::string path);
            void findAllOPACITYIndicies(GridData const *gridData);
            void findAllFeosIndicies(GridData const *gridData);
            double linearInterpolation(double lower_var, double upper_var, double corr_val, double lower_val, double upper_val);
            double biLinearInterpolation(double lowerTemp, double upperTemp, double lowerDensity, double upperDensity, double q11, double q12, double q21, 
                                    double q22, double temp, double density); 

            std::vector<double> OpacityTemperature;
            std::vector<double> OpacityDensity;
            std::vector<double> PhotonGrid;
            std::vector<double> rosslandOpacity;
            std::vector<double> planckEmissionOpacity;
            std::vector<double> planckAbsorptionOpacity;
            std::vector<double> intEnergyElectronFeosTable; 
            std::vector<double> pressureElectronFeosTable; 
            std::vector<double> intEnergyIonFeosTable; 
            std::vector<double> pressureIonFeosTable; 
            std::vector<double> ionisationFeosTable; 
            std::vector<double> FEOSTemperature;
            std::vector<double> FEOSDensity;
            std::vector<searchIndicies> searchQuantsElectron;
            std::vector<searchIndicies> searchQuantsOpacity;
            std::vector<searchIndicies> searchQuantsIon;




};
#endif