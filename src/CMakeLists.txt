add_library(mrcpp "")

add_dependencies(mrcpp mrcpp-info)

target_compile_definitions(mrcpp
  INTERFACE
    $<INSTALL_INTERFACE:USING_MRCPP>
  )

if (TARGET OpenMP::OpenMP_CXX)
  target_compile_definitions(mrcpp PUBLIC MRCPP_HAS_OMP)
  set_target_properties(mrcpp PROPERTIES MRCPP_HAS_OMP TRUE)
else()
  set_target_properties(mrcpp PROPERTIES MRCPP_HAS_OMP FALSE)
endif()

if (TARGET MPI::MPI_CXX)
  target_compile_definitions(mrcpp PUBLIC MRCPP_HAS_MPI)
  set_target_properties(mrcpp PROPERTIES MRCPP_HAS_MPI TRUE)
else()
  set_target_properties(mrcpp PROPERTIES MRCPP_HAS_MPI FALSE)
endif()

set_target_properties(mrcpp PROPERTIES EXPORT_PROPERTIES "MRCPP_HAS_OMP;MRCPP_HAS_MPI")

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
    $<BUILD_INTERFACE:$<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>>
    $<BUILD_INTERFACE:$<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>>
  )

# --- Optional FFTW3 backend (src/CMakeLists.txt) -----------------------------
option(MRCPP_USE_FFTW "Link MRCPP to FFTW3 (optional)" OFF)

if(MRCPP_USE_FFTW)
  # Try official CMake package first
  find_package(FFTW3 QUIET COMPONENTS fftw3 fftw3_threads)

  # Fallback: pkg-config (common on Homebrew/Linux)
  if(NOT TARGET FFTW3::fftw3)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(FFTW3_PC QUIET IMPORTED_TARGET GLOBAL fftw3)
      pkg_check_modules(FFTW3_THR_PC QUIET IMPORTED_TARGET GLOBAL fftw3_threads)
    endif()
  endif()

  # Link to the mrcpp target
  if(TARGET FFTW3::fftw3)
    target_link_libraries(mrcpp PUBLIC FFTW3::fftw3)
    if(TARGET FFTW3::fftw3_threads)
      find_package(OpenMP QUIET)
      if(OpenMP_CXX_FOUND)
        target_link_libraries(mrcpp PUBLIC FFTW3::fftw3_threads OpenMP::OpenMP_CXX)
        target_compile_definitions(mrcpp PUBLIC MRCPP_FFTW_THREADS)
      endif()
    endif()
  elseif(FFTW3_PC_FOUND)
    target_link_libraries(mrcpp PUBLIC PkgConfig::FFTW3_PC)
    if(FFTW3_THR_PC_FOUND)
      find_package(OpenMP QUIET)
      if(OpenMP_CXX_FOUND)
        target_link_libraries(mrcpp PUBLIC PkgConfig::FFTW3_THR_PC OpenMP::OpenMP_CXX)
        target_compile_definitions(mrcpp PUBLIC MRCPP_FFTW_THREADS)
      endif()
    endif()
  else()
    message(FATAL_ERROR "MRCPP_USE_FFTW=ON but FFTW3 was not found (CMake and pkg-config).")
  endif()

  # Expose a define so code can #ifdef it
  target_compile_definitions(mrcpp PUBLIC MRCPP_HAVE_FFTW)
endif()

add_subdirectory(core)
add_subdirectory(functions)
add_subdirectory(operators)
add_subdirectory(treebuilders)
add_subdirectory(trees)
add_subdirectory(utils)

set_target_properties(mrcpp
  PROPERTIES
    VERSION ${MRCPP_VERSION_MAJOR}
    SOVERSION ${MRCPP_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "mrcpp"
    EXPORT_NAME "mrcpp"
    PUBLIC_HEADER "${mrcpp_headers}"
  )

if(ENABLE_TESTS AND APPLE)
  set_target_properties(mrcpp
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

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
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRCPPConfigVersion.cmake
  VERSION ${MRCPP_VERSION}
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
