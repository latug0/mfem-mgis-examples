set(MFEMMGIS_EX1 "${MFEMMGIS_DIR}/../examples/ex1")
set(REF_DIR "${CMAKE_SOURCE_DIR}/ex2")
file(COPY
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${REF_DIR}/ssna303.msh
  DESTINATION ${CMAKE_BINARY_DIR}/ex2)

install(FILES
  ${REF_DIR}/CMakeLists.txt
  ${MFEMMGIS_EX1}/Plasticity.mfront
  ${REF_DIR}/ssna303.cxx
  ${REF_DIR}/ssna303.msh
  DESTINATION share/mfem-mgis-examples/ex2)


