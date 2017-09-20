if(NOT DEFINED ${MW_FILTER_DIR})
    set(MW_FILTER_DIR "${PROJECT_SOURCE_DIR}/src/mrcpp/mwfilters" CACHE STRING "Path to MW filters and cross-correlation coefs")
endif()

option(MRCPP_INSTALL_FILTERS "Install MW filters" OFF)

if(MRCPP_INSTALL_FILTERS)
    install(DIRECTORY ${MW_FILTER_DIR}
        DESTINATION share
        )
    set(MW_FILTER_DIR ${CMAKE_INSTALL_PREFIX}/share/MultiwaveletFilters)
endif()
