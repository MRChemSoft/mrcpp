target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/scaling_factor.cpp
    )

add_Catch_test(NAME scaling_factor         LABELS scaling_factor)
