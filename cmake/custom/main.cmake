include(GNUInstallDirs)

configure_file (
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )

include_directories(${PROJECT_BINARY_DIR})

include_directories(${PROJECT_BINARY_DIR}/external/include)
link_directories(${PROJECT_BINARY_DIR}/external/lib)

include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
