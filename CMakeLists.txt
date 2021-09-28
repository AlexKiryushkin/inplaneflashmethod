
cmake_minimum_required(VERSION 3.10)
project(ipfm)

list(APPEND MySources "./src/calculate_unknown_diffusivity.cpp")

add_library( ipfm_lib SHARED ${MySources})

target_compile_definitions(ipfm_lib PUBLIC ipfm_lib)
target_compile_definitions(ipfm_lib PRIVATE IPFM_EXPORTS)
target_compile_features(ipfm_lib PRIVATE cxx_std_17)
target_include_directories(ipfm_lib PRIVATE "./include")
target_include_directories(ipfm_lib PRIVATE "./deps/include")

add_executable( ipfm "./src/main.cpp"
                     "./src/read_parameters.cpp" )
target_compile_features(ipfm PRIVATE cxx_std_17)
target_link_libraries( ipfm PUBLIC ipfm_lib )
target_include_directories(ipfm PRIVATE "./include")
target_include_directories(ipfm PRIVATE "./deps/include")

add_custom_command(TARGET ipfm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/parameters.txt"
                                     "${CMAKE_BINARY_DIR}/$<CONFIG>/parameters.txt"
)
