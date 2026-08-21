cmake_minimum_required(VERSION 3.10)
project(ProjetSimuThermoMeca VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MFEMMGIS REQUIRED)

if(CMAKE_CONFIGURATION_TYPES)
if (NOT TARGET check)
add_custom_target(check COMMAND
${CMAKE_CTEST_COMMAND} -T test -C $)
endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
if (NOT TARGET check)
add_custom_target(check COMMAND
${CMAKE_CTEST_COMMAND} -T test )
endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)


get_filename_component(MFRONT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MFront" ABSOLUTE)

file(GLOB MFRONT_LOIS "${MFRONT_INCLUDE_DIR}/*.mfront")

set(MFRONT_OUTPUTS
    "${CMAKE_CURRENT_BINARY_DIR}/src/libU3SI2-generic.so"
    "${CMAKE_CURRENT_BINARY_DIR}/src/libALFENI-generic.so"
)

# Compiling MFront laws
add_custom_command(
    OUTPUT ${MFRONT_OUTPUTS}
    COMMAND mfront --obuild --interface=generic --include=${MFRONT_INCLUDE_DIR} ${MFRONT_LOIS}
    DEPENDS ${MFRONT_LOIS}
    COMMENT "Generating laws for thermomechanical simulation."
    VERBATIM
)

add_executable(Thermomechanical
    Thermomechanical.cxx 
    src/BoundaryConditions.cxx 
    src/RobinBC.cxx
)

target_include_directories(Thermomechanical PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")

add_custom_target(
    generate-libBehaviour ALL
    DEPENDS ${MFRONT_OUTPUTS}
)

target_link_libraries(Thermomechanical PRIVATE
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
	mfem-mgis::MFEMMGIS
    tbb
)

add_dependencies(check Thermomechanical generate-libBehaviour)