set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${MFEMMGIS_EX1}/cube.mesh
  ${MFEMMGIS_EX1}/UniaxialTensileTest.cxx
  ${MFEMMGIS_EX1}/UnitTestingUtilities.hxx
  DESTINATION ${CMAKE_BINARY_DIR}/ex1)

#set(MFEMMGIS_TESTDIR "${MFEMMGIS_DIR}/../examples/ex1")
#file(COPY ${MFEMMGIS_TESTDIR}/Plasticity.mfront
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#file(COPY ${MFEMMGIS_TESTDIR}/cube.mesh
#  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#install(FILES
#  ex1/CMakeLists.txt
#  ex1/Makefile
#  DESTINATION share/mfem-mgis/examples/ex1)
#
