target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi.cpp
)

list(APPEND mwparallel_h
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi.h
    ${CMAKE_CURRENT_SOURCE_DIR}/omp.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list
foreach(_header ${mwparallel_h})
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/mwparallel
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/mwparallel
    )
endforeach() 

