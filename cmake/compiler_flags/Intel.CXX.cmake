if(NOT DEFINED ENV{CXXFLAGS})
    # Intel Classic (icc/icpc) only — IntelLLVM (icx/icpx) uses Clang flags
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz -fp-speculation=fast -fp-model fast -Wno-unknown-pragmas")
        set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -debug -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG")
    endif()
endif()
