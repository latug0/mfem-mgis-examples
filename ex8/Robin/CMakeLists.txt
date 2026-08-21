cmake_minimum_required(VERSION 3.13)
project(Test_robin_thermique)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require MFEM-MGIS
find_package(MFEMMGIS REQUIRED)

# Compiling MFront laws
add_custom_command(
    OUTPUT 
        "${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.so"
    COMMAND mfront --obuild --interface=generic ../StationaryLinearHeatTransfer.mfront
    DEPENDS StationaryLinearHeatTransfer.mfront
    COMMENT "Generating laws for Robin boundary conditions test."
    VERBATIM
)

add_custom_target(
    behaviour ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.so"
)

# Define the executable to build for this example
add_executable(robin_test robin_test.cxx ../src/RobinBC.cxx)

add_dependencies(robin_test behaviour)

target_link_libraries(robin_test mfem-mgis::MFEMMGIS tbb)
