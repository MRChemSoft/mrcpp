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

add_executable(generate_filters_sources generate_filters_sources.cpp)
target_link_libraries(generate_filters_sources
  PUBLIC
    Eigen3::Eigen
  )
set_target_properties(generate_filters_sources
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

# Prepare install tree
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mwfilters
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
  )
