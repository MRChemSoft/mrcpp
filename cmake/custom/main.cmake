# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

configure_file (
  ${PROJECT_SOURCE_DIR}/config.h.in
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/config.h
  )

install(
  FILES ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(pilot)
add_subdirectory(api)
