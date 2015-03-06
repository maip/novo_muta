cmake_minimum_required(VERSION 2.6)
project(bam_driver)
project(bin_driver)
project(indiv_driver)
project(parse_bam_header)
project(pileup_driver)
project(simulation_driver)
project(simulation_trio)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated-register")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")
set(EM sufficient_statistics.cc parameter_estimates.cc)
set(TRIO_MODEL utility.cc read_dependent_data.cc trio_model.cc)

include_directories("./")
include_directories("/test")
include_directories("${CMAKE_SOURCE_DIR}/Eigen/src")

# Comment out if usr/local/include/bamtools/utils exists.
include_directories("${CMAKE_SOURCE_DIR}/bamtools/src")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
exec_program(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBRARIES)

find_package(BamTools REQUIRED)
if(BamTools_FOUND)
  include_directories(${BamTools_INCLUDE_DIRS})

  add_executable(bam_driver ${TRIO_MODEL}
                 ${BamTools_LIBRARIES} variant_visitor.cc
                 ${EM} bam_utility.cc bam_driver.cc)
endif()

add_executable(bin_driver utility.cc bin_driver.cc)
add_executable(indiv_driver ${TRIO_MODEL} indiv_driver.cc)
add_executable(parse_bam_header utility.cc parse_bam_header.cc)
add_executable(pileup_driver ${TRIO_MODEL} ${EM} pileup_utility.cc pileup_driver.cc)
add_executable(simulation_driver ${TRIO_MODEL} simulation_model.cc simulation_driver.cc)
add_executable(simulation_trio ${TRIO_MODEL} simulation_trio.cc)

target_link_libraries(bam_driver bamtools bamtools-utils ${GSL_LIBRARIES})
target_link_libraries(indiv_driver ${GSL_LIBRARIES})
target_link_libraries(simulation_driver ${GSL_LIBRARIES})
target_link_libraries(simulation_trio ${GSL_LIBRARIES})
target_link_libraries(pileup_driver ${GSL_LIBRARIES})

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.36.0 COMPONENTS system filesystem unit_test_framework REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  project(test_utility)
  add_executable(test_utility utility.cc test/test_utility.cc)
  target_link_libraries(test_utility ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_trio_model)
  add_executable(test_trio_model ${TRIO_MODEL} test/test_trio_model.cc)
  target_link_libraries(test_trio_model ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_read_dependent_data)
  add_executable(test_read_dependent_data ${TRIO_MODEL} test/test_read_dependent_data.cc)
  target_link_libraries(test_read_dependent_data ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_simulation_model)
  add_executable(test_simulation_model ${TRIO_MODEL} simulation_model.cc test/test_simulation_model.cc)
  target_link_libraries(test_simulation_model ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_pileup_utility)
  add_executable(test_pileup_utility ${TRIO_MODEL} pileup_utility.cc test/test_pileup_utility.cc)
  target_link_libraries(test_pileup_utility ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_sufficient_statistics)
  add_executable(test_sufficient_statistics ${TRIO_MODEL} ${EM} test/test_sufficient_statistics.cc)
  target_link_libraries(test_sufficient_statistics ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_parameter_estimates)
  add_executable(test_parameter_estimates ${TRIO_MODEL} ${EM} test/test_parameter_estimates.cc)
  target_link_libraries(test_parameter_estimates ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_em_algorithm)
  add_executable(test_em_algorithm ${TRIO_MODEL} ${EM} test/test_parameter_estimates.cc)
  target_link_libraries(test_em_algorithm ${GSL_LIBRARIES} ${Boost_LIBRARIES})
endif()

enable_testing()

add_test(test_utility test_utility)
add_test(test_trio_model test_trio_model)
add_test(test_read_dependent_data test_read_dependent_data)
add_test(test_simulation_model test_simulation_model)
add_test(test_pileup_utility test_pileup_utility)
add_test(test_sufficient_statistics test_sufficient_statistics)
add_test(test_parameter_estimates test_parameter_estimates)
add_test(test_em_algorithm test_em_algorithm)