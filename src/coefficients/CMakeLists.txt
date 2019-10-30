file(GLOB _cpp_sources
  LIST_DIRECTORIES FALSE
  *.cpp
  )

target_sources(mrcpp
  PRIVATE
    ${_cpp_sources}
  )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
  ${CMAKE_CURRENT_SOURCE_DIR}/FilterData.h
  )

# Sets install directory for all the headers in the list
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname})
foreach(_header IN LISTS ${_dirname}_h)
  configure_file(
    ${_header}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname}
    COPYONLY
    )
  install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname}
    )
endforeach()
