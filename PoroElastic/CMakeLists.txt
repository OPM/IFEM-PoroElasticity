cmake_minimum_required(VERSION 3.5)

project(PoroElastic)

include_directories(${PROJECT_SOURCE_DIR}/.. ${ELASTICITY_DIR})

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

file(GLOB PoroEl_SOURCES *.C)
file(GLOB PoroEl_HEADERS *.h)
add_library(PoroElastic STATIC ${PoroEl_SOURCES})

list(APPEND CHECK_SOURCES ${PoroEl_SOURCES})
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
