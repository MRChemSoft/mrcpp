# Export compile commands for each file to JSON
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

file(READ "${PROJECT_SOURCE_DIR}/VERSION" MRCPP_VERSION)
string(STRIP "${MRCPP_VERSION}" MRCPP_VERSION)

string(REPLACE "." ";" VERSION_LIST ${MRCPP_VERSION})
list(GET VERSION_LIST 0 MRCPP_VERSION_MAJOR)
list(GET VERSION_LIST 1 MRCPP_VERSION_MINOR)

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
                     -DMRCPP_VERSION=${MRCPP_VERSION}
                     -DMW_FILTER_SOURCE_DIR=${MW_FILTER_SOURCE_DIR}
                     -DMW_FILTER_INSTALL_DIR=${MW_FILTER_INSTALL_DIR}
                     -P ${CMAKE_CURRENT_LIST_DIR}/binary-info.cmake
  MAIN_DEPENDENCY
    ${PROJECT_SOURCE_DIR}/version.h.in
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}
  )

# rebuild version_info.h every time
add_custom_target(
  mrcpp-info
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E touch_nocreate ${PROJECT_SOURCE_DIR}/version.h.in
  DEPENDS
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  )

# See here for the reason why: https://gitlab.kitware.com/cmake/cmake/issues/18399
set_source_files_properties(${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  PROPERTIES
    GENERATED 1
  )

list(APPEND mrcpp_headers
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/config.h
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/version.h
  )

include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)

add_subdirectory(api)
add_subdirectory(src)
add_subdirectory(pilot)
