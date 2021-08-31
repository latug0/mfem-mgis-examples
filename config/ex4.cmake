set(MFEMMGIS_EXDIR "${MFEMMGIS_DIR}/../examples")
set(LOCAL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ex4")

file(COPY
  ${MFEMMGIS_EXDIR}/ex4/Elasticity.mfront
  ${MFEMMGIS_EXDIR}/env.sh
  DESTINATION ${LOCAL_SOURCE_DIR})

install(FILES
  ${MFEMMGIS_EXDIR}/ex4/Elasticity.mfront
  ${MFEMMGIS_EXDIR}/env.sh
  ${LOCAL_SOURCE_DIR}/CMakeLists.txt
  ${LOCAL_SOURCE_DIR}/Makefile
  ${LOCAL_SOURCE_DIR}/GrainsEx.cxx
  DESTINATION share/mfem-mgis-examples/ex4)


