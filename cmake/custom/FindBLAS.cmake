# - Find a BLAS library
#
# This module will first look in BLAS_ROOT before considering the default
# system pahts.
# The linker language can be defined by setting the varable BLAS_LANG
#
# This module defines:
#
#  BLAS_INCLUDE_DIRS Where to find blas.h (or equivalent)
#  BLAS_LIBRARIES Libraries to link against to use BLAS
#  BLAS_FOUND Defined if BLAS is available
#  HAVE_BLAS To be used in #ifdefs
#
# None of the above will be defined unless BLAS can be found.
#
#=============================================================================
# Copyright 2011 Jonas Juselius <jonas.juselius@uit.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

include(MathLibFunctions)
include(FindPackageMessage)

if (EXISTS $ENV{MATH_ROOT})
    if (NOT DEFINED BLAS_ROOT})
        set(BLAS_ROOT $ENV{MATH_ROOT})
    endif()
endif()

if (EXISTS $ENV{BLAS_ROOT})
    if (NOT DEFINED BLAS_ROOT})
        set(BLAS_ROOT $ENV{BLAS_ROOT})
    endif()
endif()

# BLAS and LAPACK often go together
if (NOT DEFINED BLAS_ROOT})
    if (DEFINED LAPACK_ROOT})
        set(BLAS_ROOT ${LAPACK_ROOT})
    elseif (EXISTS $ENV{LAPACK_ROOT})
        set(BLAS_ROOT $ENV{LAPACK_ROOT})
    endif()
endif()

if (BLAS_INCLUDE_DIRS AND BLAS_LIBRARIES)
  set(BLAS_FIND_QUIETLY TRUE)
endif ()

if (NOT BLAS_FIND_COMPONENTS)
    if (DEFINED BLAS_TYPE)
        set(BLAS_FIND_COMPONENTS ${BLAS_TYPE})
    elseif(ENABLE_64BIT_INTEGERS)
        set(BLAS_FIND_COMPONENTS MKL)
    else()
        set(BLAS_FIND_COMPONENTS MKL Atlas ACML default)
    endif()
endif()

function(find_blas)
    foreach (blas ${BLAS_FIND_COMPONENTS})
        if (${blas} MATCHES "MKL")
            find_mkl()
        elseif (${blas} MATCHES "Atlas")
            find_atlas()
        elseif (${blas} MATCHES "ACML")
            find_acml()
        else()
            find_generic()
        endif()
        if (BLAS_FOUND)
            break()
        endif()
    endforeach()
endfunction()

macro(find_generic)
    if (MATH_LANG STREQUAL "C")
        find_math_header(blas cblas.h)
    endif()
    find_math_libs(blas blas)
    cache_math_result(blas generic)
endmacro()

macro(find_acml)
    set(MATH_LIBRARY_PATH_SUFFIXES libso)
    set(MATH_INCLUDE_PATH_SUFFIXES)
    if (MATH_LANG STREQUAL "C")
        find_math_header(blas cblas.h)
    endif()
    find_math_libs(blas acml)
    cache_math_result(blas ACML)
endmacro()

macro(find_atlas)
    set(MATH_LIBRARY_PATH_SUFFIXES 
        atlas atlas-base atlas-base/atlas)
    set(MATH_INCLUDE_PATH_SUFFIXES atlas)
    if (MATH_LANG STREQUAL "C")
        find_math_header(blas cblas.h)
        find_math_libs(blas cblas atlas f77blas)
    else()
        find_math_libs(blas atlas f77blas blas)
    endif()
    cache_math_result(blas Atlas) 
endmacro()

macro(find_mkl)
    set(MATH_INCLUDE_PATH_SUFFIXES)
    if (MATH_LANG STREQUAL "C")
        find_math_header(blas mkl_cblas.h)
    endif()

    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        set(MATH_LIBRARY_PATH_SUFFIXES intel64 em64t)
        if(ENABLE_64BIT_INTEGERS)
            set(blas_libs mkl_core mkl_intel_ilp64 mkl_sequential
                guide pthread m)
        else()
            set(blas_libs mkl_core mkl_intel_lp64 mkl_sequential
                guide pthread m)
        endif()
    else()
        set(MATH_LIBRARY_PATH_SUFFIXES ia32 32)
        set(blas_libs mkl_core mkl_intel mkl_sequential guide pthread m)
    endif()
    find_math_libs(blas ${blas_libs})

    # Newer MKL versions don't have libguide
    if (NOT BLAS_LIBRARIES)
        if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
            set(MATH_LIBRARY_PATH_SUFFIXES intel64 em64t)
            if(ENABLE_64BIT_INTEGERS)
                set(blas_libs mkl_core mkl_intel_ilp64 mkl_sequential
                    pthread m)
            else()
                set(blas_libs mkl_core mkl_intel_lp64 mkl_sequential
                    pthread m)
            endif()
        else()
            set(MATH_LIBRARY_PATH_SUFFIXES ia32 32)
            set(blas_libs mkl_core mkl_intel mkl_sequential guide pthread m)
        endif()
        find_math_libs(blas ${blas_libs})
    endif()

    if (BLAS_LIBRARIES)
        set(BLAS_LIBRARIES -Wl,--start-group ${BLAS_LIBRARIES} -Wl,--end-group)
    endif()
    cache_math_result(blas MKL)
    unset(blas_libs)
endmacro()

find_blas()
if(BLAS_FOUND)
   find_package_message(BLAS "Found BLAS: ${BLAS_TYPE}" "[${BLAS_LIBRARIES}]")
endif()

