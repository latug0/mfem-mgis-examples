set(MFEMMGIS_EXDIR "${MFEMMGIS_DIR}/../examples")
set(LOCAL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ex2")

file(COPY
  ${MFEMMGIS_EXDIR}/ex1/Plasticity.mfront
  ${MFEMMGIS_EXDIR}/env.sh
  ${LOCAL_SOURCE_DIR}/ssna303.msh
  DESTINATION ${CMAKE_BINARY_DIR}/ex2)

install(FILES
  ${MFEMMGIS_EXDIR}/ex1/Plasticity.mfront
  ${MFEMMGIS_EXDIR}/env.sh
  ${LOCAL_SOURCE_DIR}/CMakeLists.txt
  ${LOCAL_SOURCE_DIR}/Makefile
  ${LOCAL_SOURCE_DIR}/ssna303.cxx
  ${LOCAL_SOURCE_DIR}/ssna303.msh
  DESTINATION share/mfem-mgis-examples/ex2)


