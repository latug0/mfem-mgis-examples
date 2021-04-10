set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${MFEMMGIS_EX1}/cube.mesh
  ${MFEMMGIS_EX1}/UniaxialTensileTest.cxx
  ${MFEMMGIS_EX1}/UnitTestingUtilities.hxx
  DESTINATION ${PROJECT_BINARY_DIR})

#set(MFEMMGIS_TESTDIR "${MFEMMGIS_DIR}/../examples/ex1")
#file(COPY ${MFEMMGIS_TESTDIR}/Plasticity.mfront
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#file(COPY ${MFEMMGIS_TESTDIR}/cube.mesh
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

