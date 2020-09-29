add_library(mrcpp "")

add_dependencies(mrcpp mrcpp-info)

target_compile_definitions(mrcpp
  INTERFACE
    $<INSTALL_INTERFACE:USING_MRCPP>
  PUBLIC
    $<$<TARGET_EXISTS:MPI::MPI_CXX>:HAVE_MPI>
  )

target_include_directories(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
  )

target_link_libraries(mrcpp
  PUBLIC
    Eigen3::Eigen
    $<BUILD_INTERFACE:$<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>>
    $<BUILD_INTERFACE:$<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>>
  )

add_subdirectory(core)
add_subdirectory(functions)
add_subdirectory(operators)
add_subdirectory(treebuilders)
add_subdirectory(trees)
add_subdirectory(utils)

set_target_properties(mrcpp
  PROPERTIES
    VERSION ${MRCPP_VERSION_MAJOR}
    SOVERSION ${MRCPP_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "mrcpp"
    EXPORT_NAME "mrcpp"
    PUBLIC_HEADER "${mrcpp_headers}"
  )

if(ENABLE_TESTS AND APPLE)
  set_target_properties(mrcpp
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

install(
  TARGETS
    mrcpp
  EXPORT
    "MRCPPTargets"
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )

# <<<  Export Config  >>>
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfigVersion.cmake
  VERSION ${MRCPP_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/MRCPPConfig.cmake.in
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfig.cmake
  INSTALL_DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

install(
  FILES
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfig.cmake
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfigVersion.cmake
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

install(
  EXPORT
    "MRCPPTargets"
  NAMESPACE
    "MRCPP::"
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
    )

# Alias useful for use with FetchContent
add_library(MRCPP::mrcpp ALIAS mrcpp)
