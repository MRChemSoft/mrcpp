project(trashy LANGUAGES CXX)
find_package(MRCPP CONFIG REQUIRED)

get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
message(STATUS "${_loc}")
