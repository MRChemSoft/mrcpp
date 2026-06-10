#.rst:
#
# Enables architecture-specific compiler flags.
#
# Variables used::
#
#   ENABLE_ARCH_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--arch-flags=<ARCH_FLAGS> Enable architecture-specific compiler flags [default: True]."
#   define: "'-DENABLE_ARCH_FLAGS={0}'.format(arguments['--arch-flags'])"

option(ENABLE_ARCH_FLAGS "Enable architecture-specific compiler flags" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Enforce minimum compiler versions for C++17 support
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "11.2")
    message(FATAL_ERROR "GCC >= 11.2 required, found ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14.0")
    message(FATAL_ERROR "Clang >= 14.0 required, found ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "2022.1")
    message(FATAL_ERROR "Intel oneAPI (icx) >= 2022.1 required, found ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
#   if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.0")
#     message(FATAL_ERROR "Intel Classic (icc) >= 19.0 required, found ${CMAKE_CXX_COMPILER_VERSION}")
#   endif()
endif()

if(ENABLE_ARCH_FLAGS)
  # iterate over list of flags and use the first one that is compatible with the
  # compiler in use
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(_arch_flag "-march=native")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    set(_arch_flag "-march=native")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(_arch_flag "-xHost")
  endif()
  message(STATUS "Adding architecture-specific compiler flag: ${_arch_flag}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_arch_flag}")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/GNU.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Intel.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Clang.CXX.cmake)

string(REPLACE " " ";" _cmake_cxx_flags ${CMAKE_CXX_FLAGS})
string(REPLACE " " ";" _cmake_cxx_flags_release ${CMAKE_CXX_FLAGS_RELEASE})
string(REPLACE " " ";" _cmake_cxx_flags_relwithdebinfo ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
