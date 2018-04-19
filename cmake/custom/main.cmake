configure_file (
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )

install(
    FILES ${CMAKE_BINARY_DIR}/config.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(pilot)
add_subdirectory(api)
