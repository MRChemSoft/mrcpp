target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/integral.cpp
    )

add_Catch_test(NAME integral         LABELS integral)
