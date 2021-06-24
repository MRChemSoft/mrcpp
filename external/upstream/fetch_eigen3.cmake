find_package(Eigen3 3.3 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )
if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
else()
  message(STATUS "Suitable Eigen3 could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Declare(eigen3_sources
    QUIET
    URL
      https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    PATCH_COMMAND patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/eigen-config-cmake.patch
    )

  FetchContent_GetProperties(eigen3_sources)

  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

  if(NOT eigen3_sources_POPULATED)
    FetchContent_Populate(eigen3_sources)

    add_subdirectory(
      ${eigen3_sources_SOURCE_DIR}
      ${eigen3_sources_BINARY_DIR}
      )
  endif()
  # Provide an alias, so linking to Eigen looks the same regardless if it was
  # found on the system or if it was fetched at configuration
  add_library(Eigen3::Eigen ALIAS eigen)
endif()
