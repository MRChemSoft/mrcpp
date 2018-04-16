target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/band_width.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bounding_box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/function_tree.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/gauss_quadrature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/derivative_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/clear_grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/build_grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/helmholtz_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_convolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/legendre_poly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mw_filter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/addition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiplication.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/projection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mw_tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/node_box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/node_index.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/poisson_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scaling_basis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_io.cpp
    )

add_Catch_test(NAME band_width              LABELS band_width)
add_Catch_test(NAME bounding_box            LABELS bounding_box)
add_Catch_test(NAME function_tree           LABELS function_tree)
add_Catch_test(NAME derivative              LABELS derivative)
add_Catch_test(NAME gauss_quadrature        LABELS gauss_quadrature)
add_Catch_test(NAME clear_grid              LABELS clear_grid)
add_Catch_test(NAME build_grid              LABELS build_grid)
add_Catch_test(NAME helmholtz_operator      LABELS helmholtz_operator)
add_Catch_test(NAME identity_convolution    LABELS identity_convolution)
add_Catch_test(NAME legendre_poly           LABELS legendre_poly)
add_Catch_test(NAME mw_filter               LABELS mw_filter)
add_Catch_test(NAME addition                LABELS addition)
add_Catch_test(NAME multiplication          LABELS multiplicationx)
add_Catch_test(NAME projection              LABELS projection)
add_Catch_test(NAME mw_tree                 LABELS mw_tree)
add_Catch_test(NAME node_box                LABELS node_box)
add_Catch_test(NAME node_index              LABELS node_index)
add_Catch_test(NAME poisson_operator        LABELS poisson_operator)
add_Catch_test(NAME polynomials             LABELS polynomials)
add_Catch_test(NAME scaling_basis           LABELS scaling_basis)
add_Catch_test(NAME tree_io                 LABELS tree_io)
