set(MFEMMGIS_EXDIR "${MFEMMGIS_DIR}/../examples")
set(LOCAL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ex3")

file(COPY
  ${MFEMMGIS_EXDIR}/ex3/Elasticity.mfront
  ${MFEMMGIS_EXDIR}/ex3/cube_2mat_per.mesh
  ${MFEMMGIS_EXDIR}/env.sh
  DESTINATION ${LOCAL_SOURCE_DIR})

install(FILES
  ${MFEMMGIS_EXDIR}/ex3/Elasticity.mfront
  ${MFEMMGIS_EXDIR}/ex3/cube_2mat_per.mesh
  ${MFEMMGIS_EXDIR}/env.sh
  ${LOCAL_SOURCE_DIR}/CMakeLists.txt
  ${LOCAL_SOURCE_DIR}/Makefile
  ${LOCAL_SOURCE_DIR}/InclusionsEx.cxx
  DESTINATION share/mfem-mgis-examples/ex3)


