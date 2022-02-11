# file(GLOB HEADER_LIST CONFIGURE_DEPEND$S "${HyKiCT_SOURCE_DIR}/include/HyKICT/*.hpp")
# file(GLOB SRC_LIST CONFIGURE_DEPENDS "${HyKiCT_SOURCE_DIR}${HyKiCT_SOURCE_DIR}/src/*.cpp")

add_library(HyKiCT_lib 
            ${HyKiCT_SOURCE_DIR}/src/Data.cpp
            ${HyKiCT_SOURCE_DIR}/src/Switches.cpp
            ${HyKiCT_SOURCE_DIR}/src/FixedData.cpp
            ${HyKiCT_SOURCE_DIR}/src/GridData.cpp
            ${HyKiCT_SOURCE_DIR}/src/EoS.cpp
            ${HyKiCT_SOURCE_DIR}/src/IdealGasEoS.cpp
            ${HyKiCT_SOURCE_DIR}/src/Ionization.cpp
            ${HyKiCT_SOURCE_DIR}/src/LoadInEoS.cpp
            ${HyKiCT_SOURCE_DIR}/src/FluidDynamics.cpp
            ${HyKiCT_SOURCE_DIR}/src/MultiRadTrans.cpp
            ${HyKiCT_SOURCE_DIR}/src/IntSource.cpp
            ${HyKiCT_SOURCE_DIR}/src/ExtSource.cpp
            ${HyKiCT_SOURCE_DIR}/src/PlasmaParameters.cpp
            ${HyKiCT_SOURCE_DIR}/src/TimeStep.cpp
            ${HyKiCT_SOURCE_DIR}/src/IO.cpp
            ${HyKiCT_SOURCE_DIR}/src/Init.cpp
            ${HyKiCT_SOURCE_DIR}/src/LoadInTables.cpp
            ${HyKiCT_SOURCE_DIR}/src/SNBSource.cpp
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/Data.hpp
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/FastMaths.hpp
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/Switches.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/FixedData.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/GridData.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/EoS.hpp
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/IdealGasEoS.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/Ionization.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/LoadInEoS.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/FluidDynamics.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/RadTrans.hpp
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/MultiRadTrans.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/IntSource.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/ExtSource.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/PlasmaParameters.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/TimeStep.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/IO.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/Init.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/LoadInTables.h
            ${HyKiCT_SOURCE_DIR}/include/HyKiCT/SNBSource.h
            )
            
add_library(HYKICT::HyKiCT_lib ALIAS HyKiCT_lib)
target_compile_features(HyKiCT_lib PRIVATE cxx_std_17)
target_compile_options(HyKiCT_lib PRIVATE -fno-math-errno)
# target_link_libraries(HyKiCT_lib PUBLIC ${PETSC_LIBRARIES} ${YAML_CPP_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(HyKiCT_lib PUBLIC  ${YAML_CPP_LIBRARIES})