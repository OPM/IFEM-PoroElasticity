cmake_minimum_required(VERSION 3.22)

project(PoroElastic LANGUAGES C CXX)

ifem_add_library(
  NAME
    PoroElastic
  ALIAS
    IFEM::PoroElasticity
  ARCHIVE_NAME
    PoroElasticity
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
    IFEM::Elasticity
)
