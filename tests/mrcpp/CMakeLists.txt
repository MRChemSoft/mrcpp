add_library(mrcpp-tests OBJECT band_width.cpp
                               bounding_box.cpp
                               function_tree.cpp
#                               gauss_quadrature.cpp
                               derivative_operator.cpp
                               clear_grid.cpp
                               build_grid.cpp
                               helmholtz_operator.cpp
                               identity_convolution.cpp
                               legendre_poly.cpp
                               mw_filter.cpp
                               addition.cpp
                               multiplication.cpp
                               projection.cpp
                               mw_tree.cpp
                               node_box.cpp
                               node_index.cpp
                               poisson_operator.cpp
                               polynomial.cpp
                               scaling_basis.cpp
                               tree_io.cpp
                               )

add_Catch_test(band_width "band_width")
add_Catch_test(bounding_box "bounding_box")
add_Catch_test(function_tree "function_tree")
add_Catch_test(derivative "derivative")
#add_Catch_test(gauss_quadrature "gauss_quadrature")
add_Catch_test(clear_grid "clear_grid")
add_Catch_test(build_grid "build_grid")
add_Catch_test(helmholtz_operator "helmholtz_operator")
add_Catch_test(identity_convolution "identity_convolution")
add_Catch_test(legendre_poly "legendre_poly")
add_Catch_test(mw_filter "mw_filter")
add_Catch_test(addition "addition")
add_Catch_test(multiplication "multiplication")
add_Catch_test(projection "projection")
add_Catch_test(mw_tree "mw_tree")
add_Catch_test(node_box "node_box")
add_Catch_test(node_index "node_index")
add_Catch_test(poisson_operator "poisson_operator")
add_Catch_test(polynomials "polynomials")
add_Catch_test(scaling_basis "scaling_basis")
add_Catch_test(tree_io "tree_io")
