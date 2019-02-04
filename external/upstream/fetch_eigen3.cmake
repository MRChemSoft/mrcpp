include(FetchContent)

find_package(Eigen3 3.3 CONFIG QUIET)
if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
else()
  message(STATUS "Suitable Eigen3 could not be located. Fetching and building!")
  FetchContent_Declare(eigen3_sources
    GIT_REPOSITORY
      https://github.com/eigenteam/eigen-git-mirror
    GIT_TAG
      3.3.7
  )

  FetchContent_GetProperties(eigen3_sources)

  if(NOT eigen3_sources_POPULATED)
    FetchContent_Populate(eigen3_sources)

    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

    add_subdirectory(
      ${eigen3_sources_SOURCE_DIR}
      ${eigen3_sources_BINARY_DIR}
      )
  endif()
  # Provide an alias, so linking to Eigen looks the same regardless if it was
  # found on the system or if it was fetched at configuration
  add_library(Eigen3::Eigen ALIAS eigen)
endif()
