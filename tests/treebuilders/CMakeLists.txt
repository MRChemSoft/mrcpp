target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/clear_grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/build_grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/addition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiplication.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/projection.cpp
    )

add_Catch_test(NAME clear_grid              LABELS clear_grid)
add_Catch_test(NAME build_grid              LABELS build_grid)
add_Catch_test(NAME addition                LABELS addition)
add_Catch_test(NAME multiplication          LABELS multiplication)
add_Catch_test(NAME square                  LABELS square)
add_Catch_test(NAME projection              LABELS projection)
