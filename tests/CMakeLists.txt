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
