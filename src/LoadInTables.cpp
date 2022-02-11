#include "HyKiCT/LoadInTables.h"

LoadInTables::LoadInTables(int startIndex, int nx)
{
   mNx = nx;
   mStartIndex = startIndex; 
}
double LoadInTables::biLinearInterpolation(double lowerTemp, double upperTemp, double lowerDensity, double upperDensity, double q11, double q12, double q21, double q22, double temp, double density)
{
    double x2x1 = upperTemp - lowerTemp; //x2 - x1;
    double y2y1 = upperDensity - lowerDensity; //y2 - y1;
    double x2x = upperTemp - temp; //x2 - x;
    double y2y = upperDensity - density; //y2 - y;
    double yy1 = density - lowerDensity; //y - y1;
    double xx1 = temp - lowerTemp; //x - x1;
    double interpolated_value = 1.0 / (x2x1 * y2y1) * ( q11 * x2x * y2y +  //q11
                                                        q21 * xx1 * y2y + //q21
                                                        q12 * x2x * yy1 + //q12
                                                        q22 * xx1 * yy1   //q22
                                                        );
    
    return interpolated_value;
}
double LoadInTables::linearInterpolation(double lower_var, double upper_var, double corr_val, double lower_val, double upper_val)
{
    double weight = (corr_val - lower_val) / 
                (upper_val - lower_val);
    double interpolated_value = lower_var * (1 - weight) + upper_var * weight; 
    return interpolated_value;

}
double LoadInTables::getGradientQuantity(double lowerTemp, double upperTemp, double lowerDensity, double upperDensity, double q11, double q12, double q21, double q22, double density)
{
    double lower_rho_specific = (q12 - q11) / (upperTemp - lowerTemp);
    double upper_rho_specific = (q22 - q21) / (upperTemp - lowerTemp);
    double weight = (density - lowerDensity) / 
                (upperDensity - lowerDensity);
    double specific_heat = lower_rho_specific * (1 - weight) + upper_rho_specific * weight;
    return specific_heat;

}
void LoadInTables::mFEOSFormatFileLoader(std::string path, std::vector<double> &mVector, bool loadAll, char delimiter)
{
    //Purpose : Load in generated Feos tables. Format is row of Temperature, row of Density and then whatever parameter after.
    //Args: Path  = Path to feos table
    //      mVector = a Vector to store the content of the table 
    //      loadAll = Bool to also load in the temperature and density to its respective private member variable. This is assumed to be constant which is usually is. 
    //NOTE: mVector HAS to correspond to the appropriate FEOS table i.e. loading in Pressure E mVector = mvectorFEOSPressureE something like this. 
    //      Also, loadAll only has to be used ONCE and the values of temperature and density is stored in FEOSTemperature and FEOSDensity.

    std::ifstream fin(path);
    std::cout << path << std::endl;
    int i = 0;
    int counter;
    if (fin.is_open())
    {
        std::string strInput;
        std::string token;

        while(getline(fin,strInput))
        {
            counter = 0;
            std::stringstream iss;
            iss << strInput;

            while(getline(iss,token,delimiter))
            {
                if((i == 0) && (loadAll))
                {
                    FEOSTemperature.push_back(stod(token));
                }
    
                else if((i == 1) && (loadAll))
                {
                    FEOSDensity.push_back(stod(token));
                }
                else if((i > 1))
                {
                    mVector.push_back(stod(token));
                }
                counter++;
            }
            i++;
        }
    }
    else
    {
        std::cout << "FILE CANNOT BE LOADED. EXITING" << std::endl;
        exit(1);
    }
    
}

void LoadInTables::mTextFileLoader(std::string path, std::vector<double> &mVector, char delimiter)
{
    std::ifstream fin(path);
    std::cout << path << std::endl;
    if (fin.is_open())
    {
        std::string strInput;
        std::string token;

        while(getline(fin,strInput))
        {
            std::stringstream iss;
            iss << strInput;

            while(getline(iss,token,delimiter))
            {
                mVector.push_back(stod(token));
            }
        }
    }
    else
    {
        std::cout << "FILE CANNOT BE LOADED. EXITING" << std::endl;
        exit(1);
    }
    
}
void LoadInTables::converTemperatureToeV(double ev_to_conversion)
{
    std::transform(FEOSTemperature.begin(), FEOSTemperature.end(), FEOSTemperature.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, ev_to_conversion));
}

LoadInTables::searchIndicies LoadInTables::mIndexSearch(std::vector<double> &temperature, std::vector<double> &density, double search_temperature, double search_density)
{
    //Purpose: Searches for the indicies that correspond to a search temperature and density. Note that two values are searched for in both simply because the 
    //          the real value does not exist and thus a bilienar interpolation has to be done which requires two neighbourhing values for the search temperature/density.
    //Args:     temperature = A vector of the temperature to be searched. For all our purpose this HAS to be FEOSTemperature.
    //          density = A vector of density to be searched. For all our purpose this HAS to be FEOSDensity.
    //          search_temperature = temperature that is beign looked for.
    //          search_density = density that is beings searched for.
    //Return:   a Vector of the indicies has following format. Temperature_index_lower, temperature_index_upper, density_index_lower, density_index_upper 
    // Creates iterator search for lower and upper bounds
    //Note: There is some redundancy in this function and that is the requirment to pass the vector that is being search on. As this is always
    //      FEOS____ in these case, this is something that can be changed.


    auto temperature_upper_bound_search = std::lower_bound(temperature.begin(), temperature.end(), search_temperature );
    auto temperature_lower_bound_search = temperature_upper_bound_search - 1;
    //auto temperature_upper_bound_search = std::upper_bound(temperature.begin(), temperature.end(), search_temperature);
    //extracts the indicies 
    int temperature_lower_index = temperature_lower_bound_search - temperature.begin();
    int temperature_upper_index = temperature_upper_bound_search - temperature.begin();
    if(temperature_upper_index == temperature.size())
    {
        std::cerr << "Temperature exceed range, setting values to maximum" << std::endl;
        temperature_upper_index = temperature.size() - 1;
        temperature_lower_index = temperature_upper_index;
    }
    if(temperature_lower_index < 0)
    {
        temperature_lower_index = 0;
        temperature_upper_index = 1;
    }

    // Creates iterator search for lower and upper bounds
    auto density_upper_bound_search = std::lower_bound(density.begin(), density.end(), search_density ); 
    auto density_lower_bound_search = density_upper_bound_search - 1;

    //auto density_upper_bound_search = std::upper_bound(density.begin(), density.end(), search_density);
    //extracts the indicies 
    int density_lower_index = density_lower_bound_search - density.begin();
    int density_upper_index = density_upper_bound_search - density.begin();
    if(density_upper_index == density.size())
    {
        std::cerr << " density exceed range, setting values to maximum "<< search_density << std::endl;
        density_upper_index = density.size() - 1;
        density_lower_index = density_upper_index;
    }
    if(density_lower_index < 0)
    {
        density_lower_index = 0;
        density_upper_index = 1;
    }
    
    searchIndicies search; 
    search.te_index = std::make_pair(temperature_lower_index, temperature_upper_index);
    search.rho_index = std::make_pair(density_lower_index, density_upper_index);
    search.te_values = std::make_pair(temperature[temperature_lower_index], temperature[temperature_upper_index]);
    search.rho_values = std::make_pair(density[density_lower_index], density[density_upper_index]);
    search.combined_index[0] = density_lower_index * mNcols + temperature_lower_index;
    search.combined_index[1] = density_lower_index * mNcols + temperature_upper_index;
    search.combined_index[2] = density_upper_index * mNcols + temperature_lower_index;
    search.combined_index[3] = density_upper_index * mNcols + temperature_upper_index;
    // std::vector<int> return_indicies = {temperature_lower_index, temperature_upper_index, density_lower_index, density_upper_index};

    return(search);
}

void LoadInTables::findAllFeosIndicies(GridData const *gridData)
{
    //Purpose: Finds all the indicies for current temperature and density profiles.
    //Args: None
    //Returns: None
    //Note: This makes full use of member variables besides the passed class DataHandler.
    searchQuantsElectron.resize(mNx);
    searchQuantsIon.resize(mNx);
    for(int i = mStartIndex; i < mNx; i++)
    {
        //Temperatures being converted to eV as the FEOS uses eV
        auto search_electron_temperature = gridData->TemperatureE[i];
        auto search_ion_temperature = gridData->TemperatureI[i];
        auto search_density = gridData->Density[i];
        searchIndicies electron_info, ion_info;
        electron_info = mIndexSearch(FEOSTemperature, FEOSDensity, search_electron_temperature, search_density);
        ion_info = mIndexSearch(FEOSTemperature, FEOSDensity, search_ion_temperature, search_density);
        searchQuantsElectron[i] = electron_info; //.insert(std::end(mSearchIndiciesElectron), std::begin(electron_indicies), std::end(electron_indicies));
        searchQuantsIon[i] = ion_info; //.insert(std::end(mSearchIndiciesIon), std::begin(ion_indicies), std::end(ion_indicies));
    }
}
void LoadInTables::findAllOPACITYIndicies(GridData const *gridData)
{
    //Purpose: Finds all the indicies for current temperature and density profiles.
    //Args: None
    //Returns: None
    //Note: This makes full use of member variables besides the passed class DataHandler.
    searchQuantsOpacity.resize(mNx);
    for(int i = mStartIndex; i < mNx; i++)
    {
        //Temperatures being converted to eV as the FEOS uses eV
        auto search_electron_temperature = gridData->TemperatureE[i];
        auto search_density = gridData->Density[i];
        searchIndicies opacity_info;
        opacity_info = mIndexSearch(OpacityTemperature, OpacityDensity, search_electron_temperature, search_density);
        searchQuantsOpacity[i] = opacity_info; //.insert(std::end(mSearchIndiciesElectron), std::begin(electron_indicies), std::end(electron_indicies));
    }
}

void LoadInTables::feosEoSLoadIn(std::string path)
{
    //Purpose: Loads in all FEOS tables given a path.
    //Args: Path = Base path for all feos tables.
    //Returns: None
    //Note: Loads into private member variables.

    std::cout << "LOADING FEOS TABLES" << std::endl;
    std::string electron_int_e_path = path + "_e1.txt";
    std::cout << electron_int_e_path << "\n";
    mFEOSFormatFileLoader(electron_int_e_path, intEnergyElectronFeosTable, true);
    mNcols = FEOSTemperature.size();
    std::string ion_int_e_path = path + "_ei1.txt";
    std::cout << ion_int_e_path << "\n";
    mFEOSFormatFileLoader(ion_int_e_path, intEnergyIonFeosTable);

    std::string electron_pre_path = path + "_p1.txt";
    std::cout << electron_pre_path << "\n";
    mFEOSFormatFileLoader(electron_pre_path, pressureElectronFeosTable);
    std::string ion_pre_path = path + "_pi1.txt";
    std::cout << ion_pre_path << "\n";
    mFEOSFormatFileLoader(ion_pre_path, pressureIonFeosTable);
}
void LoadInTables::OpacityLoadIn(std::string path)
{
    //Purpose: Loads in all FEOS tables given a path.
    //Args: Path = Base path for all feos tables.
    //Returns: None
    //Note: Loads into private member variables.
    char delimiter = '\n';
    std::cout << "LOADING OPACITY TABLES" << std::endl;
    std::string Te_path = path + "Te.txt";
    mTextFileLoader(Te_path, OpacityTemperature, delimiter);
    mNcols = OpacityTemperature.size();
    std::string rho_path = path + "Rho.txt";
    mTextFileLoader(rho_path, OpacityDensity, delimiter);
    std::string planck_emission_path = path + "Planck.txt";
    mTextFileLoader(planck_emission_path, planckEmissionOpacity, delimiter);
    std::string rossland_emission_path = path + "Rossland.txt";
    mTextFileLoader(rossland_emission_path, rosslandOpacity, delimiter);
    planckAbsorptionOpacity = planckEmissionOpacity;
}
void LoadInTables::feosIonizationLoadIn(std::string path)
{
    std::string ionisation_path = path + "_z.txt";
    std::cout << ionisation_path << "\n";
    mFEOSFormatFileLoader(ionisation_path, ionisationFeosTable, true);
    mNcols = FEOSTemperature.size();
}

void LoadInTables::PhotonGridLoadIn(std::string path)
{
    std::string photonGrid = path + "Photon_grid.txt";
    std::cout << photonGrid << "\n";
    mTextFileLoader(photonGrid, PhotonGrid);
}

// void LoadInTables::setIndicies(std::vector<int> electron_indicies, std::vector<int> ion_indicies)
// {
//     mSearchIndiciesElectron = electron_indicies;
//     mSearchIndiciesIon = ion_indicies; 
// }
// std::vector<int>  LoadInTables::getElectronIndex()
// {
//     return mSearchIndiciesElectron;
// }
// std::vector<int> LoadInTables::getIonIndex()
// {
//     return mSearchIndiciesIon;
// }



