# HyKiCT_Public_release

## 1D Lagrangian Radiation-Hydrodynamics Code utilised for Kinetic Coupling Tests

## Prerequisite

- yaml-cpp v0.5 >
- Cmake v3.13 >

## BUILD

- Create a build folder (mkdir build) and cd into build.
- Do cmake .. 
- Once cmake finishes, simply type make and everything will build. 
- Once everything has compiled test with. TEST_PATH= "/HyKiCT/tests/INTEGRATED_TESTS_INITS' ./tests/test_driver (assuming you are in the build directory) 
- IF this is successfull, you have succesfully compiled the program.

## Usage 

- Update Config.yml file to your needs.
- init files required are coord.txt(m), mass.txt(Kgm^-2), temperatureE.txt(Electron Temperature in K), temperatureI.txt(Ion Temperature in K) and density.txt(kgm^-3), Material Z and Ar.
- Additional init files include: qe.txt(multipliers or load in div.q), Load in Rad Energy Density(Energy density in space for all groups units Jm^-3). 
- Non-standard EOS/Opacity load in currently the EOS requires the file to be structured comma seperated and with Te,rho, EOS ordering. The opacity requires seperate files for Te,Rho,Photon grid and rossland/planck opacities. This will be fixed sometime in the future.  
- Once this has been initialised into some folder and config has been configured simply run exe. From base directory ("~/HyKiCT") run ./build/HyKiCT -p config.yml 
- For verbose output run with -Vb flag.
