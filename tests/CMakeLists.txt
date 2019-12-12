# Prepare "Catch" library for other executables, if not already available
if(NOT TARGET Catch)
  add_library(Catch INTERFACE)
  target_include_directories(Catch INTERFACE ${PROJECT_SOURCE_DIR}/external/catch)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(mrcpp-tests unit_tests.cpp)
add_subdirectory(treebuilders)
add_subdirectory(core)
add_subdirectory(functions)
add_subdirectory(operators)
add_subdirectory(trees)
add_subdirectory(scalingfactor)

target_link_libraries(mrcpp-tests
  mrcpp
  Catch
  )

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" _RPATH)

set_target_properties(mrcpp-tests
  PROPERTIES
  MACOSX_RPATH ON
  SKIP_BUILD_RPATH OFF
  BUILD_WITH_INSTALL_RPATH OFF
  INSTALL_RPATH "${_RPATH}"
  INSTALL_RPATH_USE_LINK_PATH ON
  )

install(
  TARGETS
  mrcpp-tests
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT bin
  )
