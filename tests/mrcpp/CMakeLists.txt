add_library(mrcpp-tests OBJECT band_width.cpp
                               bounding_box.cpp
                               function_tree.cpp
                               gauss_quadrature.cpp
                               grid_cleaner.cpp
                               grid_generator.cpp
                               helmholtz_operator.cpp
                               identity_convolution.cpp
                               legendre_poly.cpp
                               mw_adder.cpp
                               mw_filter.cpp
                               mw_multiplier.cpp
                               mw_projector.cpp
                               mw_tree.cpp
                               node_box.cpp
                               node_index.cpp
                               poisson_operator.cpp
                               polynomial.cpp
                               scaling_basis.cpp
                               )

add_Catch_test(band_width "band_width")
add_Catch_test(bounding_box "bounding_box")
add_Catch_test(function_tree "function_tree")
#add_Catch_test(gauss_quadrature "gauss_quadrature")
add_Catch_test(grid_cleaner "grid_cleaner")
add_Catch_test(grid_generator "grid_generator")
add_Catch_test(helmholtz_operator "helmholtz_operator")
add_Catch_test(identity_convolution "identity_convolution")
add_Catch_test(legendre_poly "legendre_poly")
add_Catch_test(mw_adder "mw_adder")
add_Catch_test(mw_filter "mw_filter")
add_Catch_test(mw_multiplier "mw_multiplier")
add_Catch_test(mw_projector "mw_projector")
add_Catch_test(mw_tree "mw_tree")
add_Catch_test(node_box "node_box")
add_Catch_test(node_index "node_index")
add_Catch_test(poisson_operator "poisson_operator")
add_Catch_test(polynomial "polynomial")
add_Catch_test(scaling_basis "scaling_basis")
