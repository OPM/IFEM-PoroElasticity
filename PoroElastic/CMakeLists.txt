project(PoroElastic)

cmake_minimum_required(VERSION 2.6)

include_directories(${PROJECT_SOURCE_DIR}/..
                    ${PROJECT_SOURCE_DIR}/../../Elasticity)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

set(PoroEl_SOURCES PoroElasticity.C
                   PoroMaterial.C
                   PoroMats.C
                   PoroSolutions.C
                   ThermoPoroElasticity.C)

add_library(PoroElastic STATIC ${PoroEl_SOURCES})
