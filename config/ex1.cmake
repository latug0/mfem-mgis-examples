set(MFEMMGIS_EXDIR "${MFEMMGIS_DIR}/../examples")
set(LOCAL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ex1")

file(COPY
  ${MFEMMGIS_EXDIR}/ex1/Plasticity.mfront
  ${MFEMMGIS_EXDIR}/ex1/cube.mesh
  ${MFEMMGIS_EXDIR}/ex1/UniaxialTensileTest.cxx
  ${MFEMMGIS_EXDIR}/ex1/UnitTestingUtilities.hxx
  ${MFEMMGIS_EXDIR}/env.sh
  DESTINATION ${CMAKE_BINARY_DIR}/ex1)

install(FILES
  ${MFEMMGIS_EXDIR}/ex1/Plasticity.mfront
  ${MFEMMGIS_EXDIR}/ex1/cube.mesh
  ${MFEMMGIS_EXDIR}/ex1/UniaxialTensileTest.cxx
  ${MFEMMGIS_EXDIR}/ex1/UnitTestingUtilities.hxx
  ${MFEMMGIS_EXDIR}/env.sh
  ${LOCAL_SOURCE_DIR}/CMakeLists.txt
  ${LOCAL_SOURCE_DIR}/Makefile
  DESTINATION share/mfem-mgis-examples/ex1)

