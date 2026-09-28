# Copy the law. Unlike file(COPY), configure_file updates the copy when the law
# is modified.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/StationaryLinearHeatTransfer.mfront
               ${CMAKE_CURRENT_BINARY_DIR}/StationaryLinearHeatTransfer.mfront
               COPYONLY)

# Compiling the MFront law
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.so"
    COMMAND mfront --obuild --interface=generic
            ${CMAKE_CURRENT_BINARY_DIR}/StationaryLinearHeatTransfer.mfront
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/StationaryLinearHeatTransfer.mfront
    COMMENT "Generating the law of the Robin boundary condition test.")

add_custom_target(generate-libBehaviour-robin ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.so")

add_executable(RobinTest robin_test.cxx ../src/RobinBC.cxx)
target_link_libraries(RobinTest PRIVATE MPI::MPI_CXX mfem-mgis::MFEMMGIS)

add_dependencies(check RobinTest generate-libBehaviour-robin)

add_test(NAME RobinTest
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        $<TARGET_FILE:RobinTest>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(RobinTest PROPERTIES PROCESSORS 2)
