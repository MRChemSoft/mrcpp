target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/derivative_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/helmholtz_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_convolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/poisson_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/schrodinger_evolution_operator.cpp
    )

add_Catch_test(NAME derivative_operator                     LABELS derivative_operator)
add_Catch_test(NAME helmholtz_operator                      LABELS helmholtz_operator)
add_Catch_test(NAME identity_convolution                    LABELS identity_convolution)
add_Catch_test(NAME poisson_operator                        LABELS poisson_operator)
add_Catch_test(NAME schrodinger_evolution_operator          LABELS schrodinger_evolution_operator)
