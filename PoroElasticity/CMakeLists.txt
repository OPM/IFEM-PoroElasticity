cmake_minimum_required(VERSION 3.22)

project(PoroElastic LANGUAGES C CXX)

ifem_add_library(
  NAME
    PoroElastic
  ALIAS
    IFEM::PoroElasticity
  SOURCES
    PoroElasticity.C
    PoroMaterial.C
    PoroMats.C
    PoroSolutions.C
    SIMPoroElasticity.C
  HEADERS
    PoroElasticity.h
    PoroMaterial.h
    PoroSolutions.h
    SIMPoroElasticity.h
  LIBRARIES
    Elasticity
)
set_target_properties(PoroElastic
  PROPERTIES
  ARCHIVE_OUTPUT_NAME
    PoroElasticity
)
