project(PoroElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

include_directories(${IFEM_INCLUDES} ../Common
                    ${PROJECT_SOURCE_DIR}/PoroElastic)

# Generate doxy in build tree
if(EXISTS ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
  CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in Doxyfile @ONLY)
endif()

if(EXISTS ${PROJECT_BINARY_DIR}/Doxyfile)
  add_custom_target(PoroElasticity-doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    COMMENT "Generating API documentation" VERBATIM)
  if(NOT TARGET doc)
    add_custom_target(doc)
  endif()
  add_dependencies(doc PoroElasticity-doc)
endif()

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh)

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR} ../Elasticity)

if(NOT TARGET PoroElastic)
  add_subdirectory(PoroElastic PoroElastic)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity)
  if(NOT TARGET Elasticity)
    add_subdirectory(../IFEM-Elasticity Elasticity)
  endif()
  include_directories(../IFEM-Elasticity)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/../Elasticity)
  if(NOT TARGET Elasticity)
    add_subdirectory(../Elasticity Elasticity)
  endif()
  include_directories(../Elasticity)
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(PoroElasticity main.C)

target_link_libraries(PoroElasticity PoroElastic Elasticity ${IFEM_LIBRARIES})

# Installation
install(TARGETS PoroElasticity DESTINATION bin COMPONENT bin)

if(INSTALL_DOXY)
  install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}/Apps/PoroElasticity
          COMPONENT doc
          PATTERN *.md5 EXCLUDE
          PATTERN *.map EXCLUDE)
endif(INSTALL_DOXY)

enable_testing()
include(IFEMTesting)

# Regression tests
if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else()
  ifem_add_test(Plaxis1DVerif.reg PoroElasticity)
  ifem_add_test(Plaxis1DVerif-nonmixed.reg PoroElasticity)
  ifem_add_test(OGSBenchmark1D.reg PoroElasticity)
  ifem_add_test(OGSBenchmark1D-nonmixed.reg PoroElasticity)
  if(LRSpline_FOUND)
    ifem_add_test(Plaxis1DVerif-LR.reg PoroElasticity)
    ifem_add_test(OGSBenchmark1D-LR.reg PoroElasticity)
  endif()
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
