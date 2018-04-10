if(NOT DEFINED ${MW_FILTER_DIR})
    set(MW_FILTER_DIR "${CMAKE_SOURCE_DIR}/share/mwfilters" CACHE STRING "Path to MW filters and cross-correlation coefs")
endif()

option(MRCPP_INSTALL_FILTERS "Install MW filters" ON)

if(MRCPP_INSTALL_FILTERS)
    install(DIRECTORY ${MW_FILTER_DIR}
        DESTINATION share/${PROJECT_NAME}/mwfilters
        )
    set(MW_FILTER_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/mwfilters)
endif()
