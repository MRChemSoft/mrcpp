option(ENABLE_TESTS "Enable test suite" ON)

if(ENABLE_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(tests)
endif()
