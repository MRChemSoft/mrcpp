set(MW_FILTER_SOURCE_DIR "${PROJECT_SOURCE_DIR}/share/mwfilters" CACHE STRING "Path to MW filters and cross-correlation coefs")
set(MW_FILTER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/mwfilters)
install(
  DIRECTORY share/mwfilters
  DESTINATION share/${PROJECT_NAME}
  )
