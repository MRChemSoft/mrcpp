add_library(mrcpp-objlib OBJECT "")
# FIXME Fix visibility of symbols
set_target_properties(mrcpp-objlib
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    #CXX_VISIBILITY_PRESET hidden
    #VISIBILITY_INLINES_HIDDEN 1
  )

# TODO
#target_compile_options(mrcpp-objlib )
#target_compile_definitions(mrcpp-objlib
#  PUBLIC
#    mrcpp_EXPORTS
#  )
target_include_directories(mrcpp-objlib
  PRIVATE
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
  )

add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(trees)
add_subdirectory(treebuilders)
add_subdirectory(functions)
add_subdirectory(operators)

if(NOT STATIC_LIBRARY_ONLY)
  add_library(mrcpp-shared SHARED $<TARGET_OBJECTS:mrcpp-objlib>)
  #target_compile_definitions(mrcpp-shared
  #  PUBLIC
  #    mrcpp_EXPORTS # Probably not needed?
  #  INTERFACE
  #    $<INSTALL_INTERFACE:USING_mrcpp>
  #  )
  target_include_directories(mrcpp-shared
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  # FIXME Fix versioning of DSO
  set_target_properties(mrcpp-shared
    PROPERTIES
      VERSION 1 #${PROJECT_VERSION_MAJOR}
      SOVERSION 1 #${PROJECT_VERSION_MAJOR}
      MACOSX_RPATH ON
      OUTPUT_NAME "mrcpp"
      EXPORT_NAME "mrcpp"
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  add_library(mrcpp-static STATIC $<TARGET_OBJECTS:mrcpp-objlib>)
  #target_compile_definitions(mrcpp-static
  #  PUBLIC
  #    mrcpp_STATIC_DEFINE
  #  INTERFACE
  #    $<INSTALL_INTERFACE:USING_mrcpp>
  #  )
  target_include_directories(mrcpp-static
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  set_target_properties(mrcpp-static
    PROPERTIES
       OUTPUT_NAME "mrcpp"
       EXPORT_NAME "mrcpp"
    )
endif()

if(STATIC_LIBRARY_ONLY)
  add_library(mrcpp ALIAS mrcpp-static)
else()
  add_library(mrcpp ALIAS mrcpp-shared)
endif()

# TODO Installation of targets
#install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
