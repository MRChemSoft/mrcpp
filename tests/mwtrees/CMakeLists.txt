target_sources(mrcpp-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/band_width.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bounding_box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/function_tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mw_tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/node_box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/node_index.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_io.cpp
    )

add_Catch_test(NAME band_width              LABELS band_width)
add_Catch_test(NAME bounding_box            LABELS bounding_box)
add_Catch_test(NAME function_tree           LABELS function_tree)
add_Catch_test(NAME mw_tree                 LABELS mw_tree)
add_Catch_test(NAME node_box                LABELS node_box)
add_Catch_test(NAME node_index              LABELS node_index)
add_Catch_test(NAME tree_io                 LABELS tree_io)
