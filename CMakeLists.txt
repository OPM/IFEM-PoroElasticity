project(PoroElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

set(ELASTICITY_DIR ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity)
if(NOT EXISTS ${ELASTICITY_DIR})
  set(ELASTICITY_DIR ${PROJECT_SOURCE_DIR}/../Elasticity)
endif()

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR}
                    ${ELASTICITY_DIR} ${PROJECT_SOURCE_DIR}/PoroElastic)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(${ELASTICITY_DIR} Elasticity)
endif()

if(NOT TARGET PoroElastic)
  add_subdirectory(PoroElastic)
endif()

add_executable(PoroElasticity main.C)

target_link_libraries(PoroElasticity PoroElastic Elasticity ${IFEM_LIBRARIES})

list(APPEND CHECK_SOURCES main.C)

# Installation
install(TARGETS PoroElasticity DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_BINARY_DIR} \\
                      ${PROJECT_SOURCE_DIR}/PoroElastic \\
                      ${ELASTICITY_DIR}")
add_doc_target(PoroElasticity PoroElasticity)

# Testing
enable_testing()
include(IFEMTesting)

# Regression tests
if(MPI_FOUND)
  # Add parallel tests here
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Plaxis1DVerif.reg PoroElasticity)
  ifem_add_test(Plaxis1DVerif-mixed_C0.reg PoroElasticity)
  ifem_add_test(Plaxis1DVerif-nonmixed.reg PoroElasticity)
  ifem_add_test(Plaxis1DVerif-nondim.reg PoroElasticity)
  ifem_add_test(OGSBenchmark1D.reg PoroElasticity)
  ifem_add_test(OGSBenchmark1D-nonmixed.reg PoroElasticity)
  ifem_add_test(Terzhagi.reg PoroElasticity)
  ifem_add_test(Terzhagi-nondim.reg PoroElasticity)
  ifem_add_test(SoilColumn3D.reg PoroElasticity)
  ifem_add_test(SoilColumn3D-nonmixed.reg PoroElasticity)
  ifem_add_test(VolFlux.reg PoroElasticity)
  if(LRSpline_FOUND)
    ifem_add_test(Plaxis1DVerif-LR.reg PoroElasticity)
    ifem_add_test(OGSBenchmark1D-LR.reg PoroElasticity)
  endif()
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(Terzhagi-restart.reg PoroElasticity 5)
  endif()
  ifem_add_vtf_test(SoilColumn3D.vreg PoroElasticity)
  ifem_add_hdf5_test(Terzhagi.hreg PoroElasticity)
  ifem_add_hdf5_test(SoilColumn3D.hreg PoroElasticity)
  ifem_add_hdf5_test(SoilColumn3D-nonmixed.hreg PoroElasticity)
endif()
list(APPEND TEST_APPS PoroElasticity)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  PoroElasticity
                  PoroElastic Elasticity ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
