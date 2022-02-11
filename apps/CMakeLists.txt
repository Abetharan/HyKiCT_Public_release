add_executable(HyKiCT main.cpp)
target_compile_features(HyKiCT PRIVATE cxx_std_17)
# target_link_libraries(HyKiCT PRIVATE HYKICT::HyKiCT_lib ${PETSC_LIBRARIES} ${Boost_LIBRARIES})
target_compile_options(HyKiCT PRIVATE -fno-math-errno)
target_link_libraries(HyKiCT PRIVATE HYKICT::HyKiCT_lib)

add_executable(OST_HyKiCT operator_split_T.cpp)
target_compile_features(OST_HyKiCT PRIVATE cxx_std_17)
target_compile_options(OST_HyKiCT PRIVATE -fno-math-errno)
target_link_libraries(OST_HyKiCT PRIVATE HYKICT::HyKiCT_lib )

get_target_property(MAIN_CFLAGS OST_HyKiCT COMPILE_OPTIONS)
# also see: COMPILE_DEFINITIONS INCLUDE_DIRECTORIES
message("-- Target compiler flags are: ${MAIN_CFLAGS}")
add_custom_command(TARGET OST_HyKiCT POST_BUILD
COMMAND echo built with the flags: ${MAIN_CFLAGS})
# add_executable(OS_HyKiCT operator_split_main.cpp)
# target_compile_features(OS_HyKiCT PRIVATE cxx_std_17)
# target_link_libraries(OS_HyKiCT PRIVATE HYKICT::HyKiCT_lib ${PETSC_LIBRARIES})
