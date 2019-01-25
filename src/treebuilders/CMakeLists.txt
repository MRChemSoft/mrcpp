target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ABGVCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CopyAdaptor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PHCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ProjectionCalculator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TreeBuilder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/add.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/apply.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiply.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/project.cpp
  )

list(APPEND treebuilders_h
  ${CMAKE_CURRENT_SOURCE_DIR}/ABGVCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AdditionCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/BandWidthAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CopyAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DefaultCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MultiplicationCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/OperatorAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PHCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PowerCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ProjectionCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SplitAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SquareCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TreeAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TreeBuilder.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TreeCalculator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/WaveletAdaptor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/add.h
  ${CMAKE_CURRENT_SOURCE_DIR}/apply.h
  ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/multiply.h
  ${CMAKE_CURRENT_SOURCE_DIR}/project.h
  )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS treebuilders_h)
 install(
   FILES ${_header}
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/treebuilders
   )
endforeach()
