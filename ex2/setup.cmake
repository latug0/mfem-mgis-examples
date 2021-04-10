set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${PROJECT_SOURCE_DIR}/ssna303.msh
  DESTINATION ${PROJECT_BINARY_DIR})

  #set(MFEMMGIS_TESTDIR "${MFEMMGIS_DIR}/../examples/ex2")
#file(COPY
#  ${MFEMMGIS_TESTDIR}/Plasticity.mfront
#  ${MFEMMGIS_TESTDIR}/ssna303.msh
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
#
