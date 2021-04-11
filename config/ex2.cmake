set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${CMAKE_SOURCE_DIR}/ex2/ssna303.msh
  DESTINATION ${CMAKE_BINARY_DIR}/ex2)

  #set(MFEMMGIS_TESTDIR "${MFEMMGIS_DIR}/../examples/ex2")
#file(COPY
#  ${MFEMMGIS_TESTDIR}/Plasticity.mfront
#  ${MFEMMGIS_TESTDIR}/ssna303.msh
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
#
#install(FILES
#  ex2/CMakeLists.txt
#  ex2/ssna303.cxx
#  ${mfem-mgis-data_SOURCE_DIR}/mesh/ssna303.msh
#  DESTINATION share/mfem-mgis/examples/ex2)


