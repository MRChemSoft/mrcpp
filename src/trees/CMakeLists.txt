target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/BandWidth.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundingBox.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionNode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionTreeVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GenNode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HilbertPath.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWNode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MultiResolutionAnalysis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NodeBox.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorNode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedNode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialFunctionTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialOperatorTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TreeIterator.cpp
)

list(APPEND trees_h
    ${CMAKE_CURRENT_SOURCE_DIR}/BandWidth.h
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundingBox.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionNode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionTreeVector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GenNode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/HilbertIterator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/HilbertPath.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LebesgueIterator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MWNode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MWTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MultiResolutionAnalysis.h
    ${CMAKE_CURRENT_SOURCE_DIR}/NodeBox.h
    ${CMAKE_CURRENT_SOURCE_DIR}/NodeIndex.h
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorNode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedNode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialFunctionTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialOperatorTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SerialTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/TreeIterator.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header IN LISTS trees_h)
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/trees
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/trees
    )
endforeach() 
