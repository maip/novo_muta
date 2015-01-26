/**
 * @file test_simulation_model.cc
 * @author Melissa Ip
 *
 * This file tests the functions in simulation_model.h.
 */
#define BOOST_TEST_MODULE TestSimulationModel

#include <boost/test/unit_test.hpp>

#include "simulation_model.h"


BOOST_AUTO_TEST_CASE(TestSimulationModel) {
  SimulationModel sim(4, 0.001, 1e-6, 1e-6);
  BOOST_CHECK(sim.coverage() == 4);
  BOOST_CHECK(sim.population_mutation_rate() == 0.001);
  BOOST_CHECK(sim.germline_mutation_rate() == 1e-6);
  BOOST_CHECK(sim.somatic_mutation_rate() == 1e-6);
  BOOST_CHECK(sim.has_mutation() == false);
}

BOOST_AUTO_TEST_CASE(TestSeed) {
  SimulationModel sim(4, 0.001, 1e-6, 1e-6);
  sim.Seed();
  BOOST_CHECK(string(gsl_rng_name(sim.generator())).compare("mt19937") == 0);
  sim.Free();
  // BOOST_CHECK(sim.generator() == NULL);
}

// To test these functions, run simulation_driver and call the appropriate
// count program (count_bin.cc or counts_probability.cc) to observe results
// manually.
// BOOST_AUTO_TEST_CASE(TestWriteProbability) {}
// BOOST_AUTO_TEST_CASE(TestWriteMutationCounts) {}
// BOOST_AUTO_TEST_CASE(TestPrintMutationCounts) {}