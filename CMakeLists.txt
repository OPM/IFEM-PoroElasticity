cmake_minimum_required(VERSION 3.22)

project(PoroElasticity LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

ifem_add_dependency_dir(APP Elasticity TARGET Elasticity)

if(NOT TARGET PoroElastic)
  add_subdirectory(PoroElastic)
endif()

ifem_add_application(
  NAME
    PoroElasticity
  SOURCES
    main.C
  HEADERS
    SIMDynPoroElasticity.h
    SIMStatPoroElasticity.h
  LIBRARIES
    IFEMAppCommon
    PoroElastic
)

# Testing
enable_testing()

ifem_add_hdf5_test(
  TARGET
    PoroElasticity
  TEST_FILES
    Terzhagi.hreg
    SoilColumn3D.hreg
    SoilColumn3D-nonmixed.hreg
)

ifem_add_vtf_test(
  TARGET
    PoroElasticity
  TEST_FILES
    SoilColumn3D.vreg
)

ifem_add_restart_test(
  TARGET
    PoroElasticity
  TEST_FILES
    Terzhagi-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    PoroElasticity
  TEST_FILES
    Plaxis1DVerif.reg
    Plaxis1DVerif-mixed_C0.reg
    Plaxis1DVerif-nonmixed.reg
    Plaxis1DVerif-nondim.reg
    OGSBenchmark1D.reg
    OGSBenchmark1D-nonmixed.reg
    Terzhagi.reg
    Terzhagi-nondim.reg
    SoilColumn3D.reg
    SoilColumn3D-nonmixed.reg
    VolFlux.reg
)

ifem_add_regression_test(
  TARGET
    PoroElasticity
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Plaxis1DVerif-LR.reg
    OGSBenchmark1D-LR.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    PoroElasticity
  SOURCES
    Test/TestPoroMaterial.C
    Test/TestSIMPoroElasticity.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
  LIBRARIES PoroElastic
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(
  TARGET
    PoroElasticity
  EXTRA_DIRS
    ${Elasticity_DIRECTORY}
)

# Installation
include(GNUInstallDirs)
install(TARGETS PoroElasticity DESTINATION ${CMAKE_INSTALL_BINDIR})
