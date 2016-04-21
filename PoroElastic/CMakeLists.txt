project(PoroElastic)

cmake_minimum_required(VERSION 2.6)

include_directories(${PROJECT_SOURCE_DIR}/..
                    ${PROJECT_SOURCE_DIR}/../../Elasticity)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

file(GLOB PoroEl_SOURCES *.C)
file(GLOB PoroEl_HEADERS *.h)
add_library(PoroElastic STATIC ${PoroEl_SOURCES})
