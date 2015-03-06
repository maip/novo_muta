/**
 * @file test_em_algorithm.cc
 * @author Melissa Ip
 *
 * This file tests the functions in m_algorithm.cc.
 */
#define BOOST_TEST_MODULE TestEmAlgorithm

#include <boost/test/unit_test.hpp>

#include "em_algorithm.h"


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
    BOOST_CHECK(params.sequencing_error_rate() == max_e);
    BOOST_CHECK(max_e >= 0.0 && max_e <= 1.0);
    BOOST_CHECK(stats->IsLogLikelihoodIncreasing());
  }

  for (auto trio : trios) {
    TrioModel params;
    TrioVector sites;
    sites.push_back(trio);

    BOOST_REQUIRE(sites.size() == 1);
    ParameterEstimates *stats = EstimateParameters(params, sites);
    BOOST_CHECK(params.sequencing_error_rate() == max_e);
    double max_e = stats->max_e();
    BOOST_CHECK(max_e >= 0.0 && max_e <= 1.0);
    BOOST_CHECK(stats->IsLogLikelihoodIncreasing());
  }
}
