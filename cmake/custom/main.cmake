# Export compile commands for each file to JSON
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

configure_file(
  ${PROJECT_SOURCE_DIR}/config.h.in
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/config.h
  @ONLY
  )

add_custom_command(
  OUTPUT
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  COMMAND
    ${CMAKE_COMMAND} -DINPUT_DIR=${PROJECT_SOURCE_DIR}
                     -DTARGET_DIR=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
                     -DCMAKE_SYSTEM=${CMAKE_SYSTEM}
                     -DCMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}
                     -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
                     -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                     -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                     -DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                     -DVERSION_FILE=${PROJECT_SOURCE_DIR}/VERSION
                     -DMW_FILTER_SOURCE_DIR=${MW_FILTER_SOURCE_DIR}
                     -DMW_FILTER_INSTALL_DIR=${MW_FILTER_INSTALL_DIR}
                     -P ${CMAKE_CURRENT_LIST_DIR}/binary-info.cmake
  DEPENDS
    ${PROJECT_SOURCE_DIR}/version.h.in
  BYPRODUCTS
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}
  )

# rebuild version_info.h every time
add_custom_target(
  binary-info
  ALL
  COMMAND
    touch ${PROJECT_SOURCE_DIR}/version.h.in
  DEPENDS
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  )

# See here for the reason why: https://gitlab.kitware.com/cmake/cmake/issues/18399
set_source_files_properties(${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  PROPERTIES
    GENERATED 1
  )

list(APPEND mrcpp_headers ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/config.h)
list(APPEND mrcpp_headers ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h)

include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)

add_subdirectory(api)
add_subdirectory(src)
add_subdirectory(pilot)
