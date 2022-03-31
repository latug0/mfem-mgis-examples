cmake_minimum_required(VERSION 3.13)
project(Config)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require MFEM
if (NOT TARGET mfem)
  find_package(MFEM REQUIRED)
endif(NOT TARGET mfem)

# Require MGIS
if (NOT TARGET mgis::MFrontGenericInterface)
  find_package(MFrontGenericInterface REQUIRED)
endif(NOT TARGET mgis::MFrontGenericInterface)

# Require MFEM-MGIS
if (NOT TARGET mfem-mgis::MFEMMGIS)
  find_package(MFEMMGIS REQUIRED)
endif(NOT TARGET mfem-mgis::MFEMMGIS)

if (NOT MFEM_USE_MUMPS AND NOT MFEM_USE_SUITESPARSE)
  message(FATAL_ERROR "\nYou shall have MUMPS or SUITESPARSE support activated in MFEM.\nCMake will exit." )
endif()

# Specific setup for preparing the exemples
include(ex1.cmake)
include(ex2.cmake)
include(ex3.cmake)

