cmake_minimum_required(VERSION 3.11)
project(polar_pfc VERSION 0.1 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(cmake/SHTNS.cmake)

find_path(NANOFLANN_INC nanoflann.hpp HINTS ${PROJECT_SOURCE_DIR} ${NANOFLANN_INC})
list(APPEND INCLUDE_DIRS ${NANOFLANN_INC})

find_package(ALGLib REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

add_library(config INTERFACE)
target_include_directories(config INTERFACE ${PROJECT_SOURCE_DIR} ${INCLUDE_DIRS})
target_compile_features(config INTERFACE cxx_std_11)
target_compile_definitions(config INTERFACE "SHTNS_HAS_COMPRESSION=1")
target_link_libraries(config INTERFACE ${LIBRARIES} ALGLib::ALGLib OpenMP::OpenMP_CXX ZLIB::ZLIB)

add_library(utils OBJECT
  ${PROJECT_SOURCE_DIR}/io/VtkTimeseriesWriter.cpp
  ${PROJECT_SOURCE_DIR}/io/VtkWriter.cpp
  ${PROJECT_SOURCE_DIR}/utility/Filesystem.cpp
  ${PROJECT_SOURCE_DIR}/utility/Timer.cpp)
target_link_libraries(utils PRIVATE config)

# define the executables
add_executable(polar_pfc polar_pfc2.cc)
target_link_libraries(polar_pfc PRIVATE utils config)

add_executable(qtensor qtensor.cpp)
target_link_libraries(qtensor PRIVATE utils config)
