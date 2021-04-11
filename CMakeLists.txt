
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


add_subdirectory(config)
add_subdirectory(ex1)
add_subdirectory(ex2)


#include(FetchContent)
#FetchContent_Declare(mfem-mgis
#  GIT_REPOSITORY https://github.com/thelfer/mfem-mgis.git
#  GIT_TAG master
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND ""
#  )
#
#FetchContent_GetProperties(mfem-mgis-data)
#if(NOT mfem-mgis_POPULATED)
#  FetchContent_Populate(mfem-mgis)
#endif()
