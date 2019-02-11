# Export compile commands for each file to JSON
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

configure_file (
  ${PROJECT_SOURCE_DIR}/config.h.in
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/config.h
  )

list(APPEND mrcpp_headers ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/config.h)

include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)

add_subdirectory(api)
add_subdirectory(src)
add_subdirectory(pilot)
