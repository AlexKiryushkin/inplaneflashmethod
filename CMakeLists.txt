
cmake_minimum_required(VERSION 3.10)
project(ipfm)

list(APPEND MySources "calculate_unknown_diffusivity.cpp")

add_library( ipfm_lib SHARED ${MySources})

target_compile_definitions(ipfm_lib PUBLIC ipfm_lib)
target_compile_definitions(ipfm_lib PRIVATE IPFM_EXPORTS)
target_compile_features(ipfm_lib PRIVATE cxx_std_17)

add_executable( ipfm main.cpp
                     read_parameters.cpp )
target_compile_features(ipfm PRIVATE cxx_std_17)
target_link_libraries( ipfm PUBLIC ipfm_lib )