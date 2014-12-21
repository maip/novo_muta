/**
 * @file test_parameter_estimates.cc
 * @author Melissa Ip
 *
 * This file tests the functions in parameter_estimates.h.
 */
#define BOOST_TEST_MODULE TestEmAlgorithm

#include <boost/test/unit_test.hpp>

#include "em_algorithm.cc"


BOOST_AUTO_TEST_CASE(TestEstimateParameters) {
  const TrioVector trios = GetTrioVector(kNucleotideCount);
  for (auto trio : trios) {
    TrioModel params;
    TrioVector sites;
    for (int i = 0; i < 10; ++i) {
      sites.push_back(trio);
    }

    BOOST_REQUIRE(sites.size() == 10);
    ParameterEstimates *stats = EstimateParameters(params, sites);
    double max_e = stats->max_e();
    BOOST_CHECK(max_e >= 0.0 && max_e <= 1.0);
    BOOST_CHECK(stats->IsLogLikelihoodIncreasing());
  }

  for (auto trio : trios) {
    TrioModel params;
    TrioVector sites;
    sites.push_back(trio);

    BOOST_REQUIRE(sites.size() == 1);
    ParameterEstimates *stats = EstimateParameters(params, sites);
    double max_e = stats->max_e();
    BOOST_CHECK(max_e >= 0.0 && max_e <= 1.0);
    BOOST_CHECK(stats->IsLogLikelihoodIncreasing());
  }
}
