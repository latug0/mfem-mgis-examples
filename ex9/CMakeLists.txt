cmake_minimum_required(VERSION 3.13)
project(ThirdMedium)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require MFEM-MGIS
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

# Add the tools library for these files
add_library(utils_lib STATIC ThirdMediumUtils.cxx)
target_link_libraries(utils_lib mfem-mgis::MFEMMGIS)


# Extract include path directly from the imported target
get_target_property(MFEMMGIS_INC mfem-mgis::MFEMMGIS INTERFACE_INCLUDE_DIRECTORIES)

# Search directly for the file without compiling test code
find_file(ADDITIONAL_CRITERION_HEADER
    NAMES MFEMMGIS/AbstractAdditionalConvergenceCriterion.hxx
    PATHS ${MFEMMGIS_INC}
    NO_DEFAULT_PATH
)

if(ADDITIONAL_CRITERION_HEADER)
    message(STATUS ">>> [MFEM-MGIS] Additional convergence criterion: ENABLED (${ADDITIONAL_CRITERION_HEADER})")
    target_compile_definitions(utils_lib PRIVATE MFEM_MGIS_HAS_ADDITIONAL_CONVERGENCE_CRITERION)
else()
    message(STATUS ">>> [MFEM-MGIS] Additional convergence criterion: DISABLED")
endif()


# Define the executable to build for this example
add_executable(BendingTest BendingTest.cxx)
target_link_libraries(BendingTest mfem-mgis::MFEMMGIS)
target_link_libraries(BendingTest utils_lib)

add_executable(BendingTest3D BendingTest3D.cxx)
target_link_libraries(BendingTest3D mfem-mgis::MFEMMGIS)
target_link_libraries(BendingTest3D utils_lib)

add_executable(CShapeContact CShapeContact.cxx)
target_link_libraries(CShapeContact mfem-mgis::MFEMMGIS)
target_link_libraries(CShapeContact utils_lib)


# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/Elasticity.mfront
  ${PROJECT_SOURCE_DIR}/ThirdMediumContactHyperElasticBehaviour.mfront
  ${PROJECT_SOURCE_DIR}/meshes/beam_bending/msh_files/bending.msh
  ${PROJECT_SOURCE_DIR}/meshes/beam_bending/msh_files/bending_quad.msh
  ${PROJECT_SOURCE_DIR}/meshes/beam_bending/msh_files/bending3D.msh
  ${PROJECT_SOURCE_DIR}/meshes/beam_bending/msh_files/bending3D_thin.msh
  ${PROJECT_SOURCE_DIR}/meshes/cshape/msh_files/cshape.msh
  ${PROJECT_SOURCE_DIR}/meshes/cshape/msh_files/cshape_full.msh
  ${PROJECT_SOURCE_DIR}/meshes/cshape/msh_files/point_load_approx.msh
  ${PROJECT_SOURCE_DIR}/meshes/cshape/msh_files/point_load_approx_full.msh
  ${PROJECT_SOURCE_DIR}/visualize_displacements.py
  DESTINATION ${PROJECT_BINARY_DIR})

# Use mfront to generate libBehaviour.so
# add_custom_command(
#   OUTPUT "src/libBehaviour.so"
#   COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/ThirdMediumContactHyperElasticBehaviour.mfront ${PROJECT_BINARY_DIR}/Elasticity.mfront DEPENDS ${PROJECT_BINARY_DIR}/Elasticity.mfront  ${PROJECT_BINARY_DIR}/ThirdMediumContactHyperElasticBehaviour.mfront
#   COMMENT "Performing mfront code generation")

# add_custom_target(generate-libBehaviour ALL
#   DEPENDS "src/libBehaviour.so")


