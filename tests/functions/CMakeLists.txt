target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/legendre_poly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gaussians.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/slater.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/periodify_gaussians.cpp
    )

add_Catch_test(NAME legendre_poly            LABELS legendre_poly)
add_Catch_test(NAME polynomials              LABELS polynomials)
add_Catch_test(NAME gaussians                LABELS gaussians)
add_Catch_test(NAME slater                   LABELS slater)
add_Catch_test(NAME periodic_narrow_gaussian LABELS periodic_narrow_gaussian)
add_Catch_test(NAME periodic_wide_gaussian   LABELS periodic_wide_gaussian)
add_Catch_test(NAME periodic_gaussExp        LABELS periodic_gausExp)
