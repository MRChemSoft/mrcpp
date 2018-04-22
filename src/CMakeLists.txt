add_library(mrcpp-objlib OBJECT "")
# FIXME Fix visibility of symbols
set_target_properties(mrcpp-objlib
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    #CXX_VISIBILITY_PRESET hidden
    #VISIBILITY_INLINES_HIDDEN 1
  )

# TODO Set compiler options at this level, rather than globally
#target_compile_options(mrcpp-objlib )
# TODO Set compile definitions at this level, rather than globally
#target_compile_definitions(mrcpp-objlib
#  PUBLIC
#    mrcpp_EXPORTS
#  )
target_include_directories(mrcpp-objlib
  PRIVATE
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
  )

add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(trees)
add_subdirectory(treebuilders)
add_subdirectory(functions)
add_subdirectory(operators)

if(NOT STATIC_LIBRARY_ONLY)
  add_library(mrcpp-shared SHARED $<TARGET_OBJECTS:mrcpp-objlib>)
  target_compile_definitions(mrcpp-shared
    #PUBLIC
    #  MRCPP_EXPORTS # Probably not needed?
    INTERFACE
      $<INSTALL_INTERFACE:USING_MRCPP>
    )
  target_include_directories(mrcpp-shared
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/MRCPP>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  # FIXME Fix versioning of DSO
  set_target_properties(mrcpp-shared
    PROPERTIES
      VERSION 1 #${PROJECT_VERSION_MAJOR}
      SOVERSION 1 #${PROJECT_VERSION_MAJOR}
      MACOSX_RPATH ON
      OUTPUT_NAME "mrcpp"
      EXPORT_NAME "mrcpp"
    )
  install(
    TARGETS
      mrcpp-shared
    EXPORT
      "MRCPPTargets-shared"
    RUNTIME DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  add_library(mrcpp-static STATIC $<TARGET_OBJECTS:mrcpp-objlib>)
  target_compile_definitions(mrcpp-static
    #PUBLIC
    #  MRCPP_STATIC_DEFINE
    INTERFACE
      $<INSTALL_INTERFACE:USING_MRCPP>
    )
  target_include_directories(mrcpp-static
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/MRCPP>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  set_target_properties(mrcpp-static
    PROPERTIES
       OUTPUT_NAME "mrcpp"
       EXPORT_NAME "mrcpp"
    )
  install(
    TARGETS
      mrcpp-static
    EXPORT
      "MRCPPTargets-static"
    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(STATIC_LIBRARY_ONLY)
  add_library(mrcpp ALIAS mrcpp-static)
else()
  add_library(mrcpp ALIAS mrcpp-shared)
endif()

# Generate header with preprocessor macros for visibility of symbols
include(GenerateExportHeader)
generate_export_header(mrcpp
  BASE_NAME "MRCPP"
  EXPORT_MACRO_NAME "MRCPP_EXPORT"
  EXPORT_FILE_NAME "MRCPPExport.h"
  DEPRECATED_MACRO_NAME "MRCPP_DEPRECATED"
  NO_EXPORT_MACRO_NAME "MRCPP_NO_EXPORT"
  STATIC_DEFINE "MRCPP_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "MRCPP_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED)
file(
  COPY ${CMAKE_CURRENT_BINARY_DIR}/MRCPPExport.h
  DESTINATION ${PROJECT_BINARY_DIR}/include
  )
install(
  FILES ${PROJECT_BINARY_DIR}/include/MRCPPExport.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
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

if(NOT STATIC_LIBRARY_ONLY)
  install(
    EXPORT
      "MRCPPTargets-shared"
    NAMESPACE
      "MRCPP::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  install(
    EXPORT
      "MRCPPTargets-static"
    NAMESPACE
      "MRCPP::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()
