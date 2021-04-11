set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
set(REF_DIR "${CMAKE_SOURCE_DIR}/ex1")

file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${MFEMMGIS_EX1}/cube.mesh
  ${MFEMMGIS_EX1}/UniaxialTensileTest.cxx
  ${MFEMMGIS_EX1}/UnitTestingUtilities.hxx
  DESTINATION ${CMAKE_BINARY_DIR}/ex1)

install(FILES
  ${REF_DIR}/CMakeLists.txt
  ${REF_DIR}/Makefile
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${MFEMMGIS_EX1}/cube.mesh
  ${MFEMMGIS_EX1}/UniaxialTensileTest.cxx
  ${MFEMMGIS_EX1}/UnitTestingUtilities.hxx
  DESTINATION share/mfem-mgis-examples/ex1)

