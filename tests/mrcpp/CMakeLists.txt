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

add_catch_test(band_width "band_width")
add_catch_test(bounding_box "bounding_box")
add_catch_test(function_tree "function_tree")
add_catch_test(derivative "derivative")
#add_catch_test(gauss_quadrature "gauss_quadrature")
add_catch_test(clear_grid "clear_grid")
add_catch_test(build_grid "build_grid")
add_catch_test(helmholtz_operator "helmholtz_operator")
add_catch_test(identity_convolution "identity_convolution")
add_catch_test(legendre_poly "legendre_poly")
add_catch_test(mw_filter "mw_filter")
add_catch_test(addition "addition")
add_catch_test(multiplication "multiplication")
add_catch_test(projection "projection")
add_catch_test(mw_tree "mw_tree")
add_catch_test(node_box "node_box")
add_catch_test(node_index "node_index")
add_catch_test(poisson_operator "poisson_operator")
add_catch_test(polynomials "polynomials")
add_catch_test(scaling_basis "scaling_basis")
add_catch_test(tree_io "tree_io")
