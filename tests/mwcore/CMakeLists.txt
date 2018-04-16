target_sources(mrcpp-tests PRIVATE
#   ${CMAKE_CURRENT_SOURCE_DIR}/gauss_quadrature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mw_filter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scaling_basis.cpp
    )

#add_Catch_test(NAME gauss_quadrature        LABELS gauss_quadrature)
add_Catch_test(NAME mw_filter               LABELS mw_filter)
add_Catch_test(NAME scaling_basis           LABELS scaling_basis)
