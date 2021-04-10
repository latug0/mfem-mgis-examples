
cmake_minimum_required(VERSION 3.13)
project(Mfem-mgis-examples)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)


add_subdirectory(ex1)
add_subdirectory(ex2)

#install(FILES
#  ex1/CMakeLists.txt
#  ex1/Makefile
#  DESTINATION share/mfem-mgis/examples/ex1)
#
#install(FILES
#  ex2/CMakeLists.txt
#  ex2/ssna303.cxx
#  ${mfem-mgis-data_SOURCE_DIR}/mesh/ssna303.msh
#  DESTINATION share/mfem-mgis/examples/ex2)


#include(FetchContent)
#FetchContent_Declare(mfem-mgis-data
#  GIT_REPOSITORY https://github.com/latug0/mfem-mgis-data.git
#  GIT_TAG master
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND ""
#  )
#
#FetchContent_GetProperties(mfem-mgis-data)
#if(NOT mfem-mgis-data_POPULATED)
#  FetchContent_Populate(mfem-mgis-data)
#endif()
