# Relative path to MW filters and cross-correlation coefs
set(MW_FILTER_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/mwfilters")

# Prepare build tree
file(GLOB _filters
  "${CMAKE_CURRENT_SOURCE_DIR}/mwfilters/*"
)
file(
  COPY ${_filters}
  DESTINATION ${PROJECT_NAME}/mwfilters
  )

# Prepare install tree
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mwfilters
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
  )
