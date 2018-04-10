add_library(mrcpp STATIC)

add_subdirectory(mwcore)
add_subdirectory(mwutils)
add_subdirectory(mwtrees)
add_subdirectory(mwbuilders)
add_subdirectory(mwfunctions)
add_subdirectory(mwoperators)
add_subdirectory(mwparallel)

target_include_directories(mrcpp PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/mwcore
    ${CMAKE_SOURCE_DIR}/src/mwutils
    ${CMAKE_SOURCE_DIR}/src/mwtrees
    ${CMAKE_SOURCE_DIR}/src/mwbuilders
    ${CMAKE_SOURCE_DIR}/src/mwfunctions
    ${CMAKE_SOURCE_DIR}/src/mwoperators
    )

install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


list(APPEND mrcpp_h 
    constants.h
    eigen_disable_warnings.h
    macros.h
    mrcpp_declarations.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mrcpp_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

