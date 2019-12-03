add_library(mrcpp "")

# TODO Set compiler options at this level, rather than globally
#target_compile_options(mrcpp )

target_compile_definitions(mrcpp
  INTERFACE
    $<INSTALL_INTERFACE:USING_MRCPP>
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
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
  )

add_subdirectory(core)
add_subdirectory(functions)
add_subdirectory(operators)
add_subdirectory(treebuilders)
add_subdirectory(trees)
add_subdirectory(utils)

# FIXME Fix versioning of DSO
set_target_properties(mrcpp
  PROPERTIES
    VERSION 1 #${PROJECT_VERSION_MAJOR}
    SOVERSION 1 #${PROJECT_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "mrcpp"
    EXPORT_NAME "mrcpp"
    PUBLIC_HEADER "${mrcpp_headers}"
  )

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
# FIXME Fix versioning
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfigVersion.cmake
  VERSION 1.0.0 #${PROJECT_VERSION}
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
