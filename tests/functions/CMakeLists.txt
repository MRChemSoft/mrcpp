target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/legendre_poly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gaussians.cpp
    )

add_Catch_test(NAME legendre_poly           LABELS legendre_poly)
add_Catch_test(NAME polynomials             LABELS polynomials)
add_Catch_test(NAME gaussians               LABELS gaussians)
