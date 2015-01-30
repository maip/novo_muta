/**
 * @file test_multinomial_trio_model.cc
 * @author Melissa Ip
 *
 * This file tests the functions in multinomial_trio_model.h.
 */
#define BOOST_TEST_MODULE TestMultinomialTrioModel

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "multinomial_trio_model.h"


BOOST_AUTO_TEST_CASE(TestMultinomialTrioModel) {
  MultinomialTrioModel default_params;
  RowVector4d frequencies;
  frequencies << 0.25, 0.25, 0.25, 0.25;
  double germline_mutation_rate = default_params.germline_mutation_rate();
  double exp_term = exp(-4.0/3.0 * germline_mutation_rate);
  double homozygous_match = 0.25 + 0.75 * exp_term;
  double heterozygous_match = 0.25 + 0.25 * exp_term;
  double mismatch = 0.25 - 0.25 * exp_term;

  BOOST_CHECK(default_params.population_mutation_rate() == 0.001);
  BOOST_CHECK(germline_mutation_rate == 2e-8);
  BOOST_CHECK(default_params.somatic_mutation_rate() == 2e-8);
  BOOST_CHECK(default_params.sequencing_error_rate() == 0.005);
  BOOST_CHECK(default_params.dirichlet_dispersion() == 1.0);
  BOOST_CHECK(default_params.nucleotide_frequencies() == frequencies);
  BOOST_CHECK(default_params.homozygous_match() == homozygous_match);
  BOOST_CHECK(default_params.heterozygous_match() == heterozygous_match);
  BOOST_CHECK(default_params.mismatch() == mismatch);
  BOOST_CHECK_CLOSE(default_params.population_priors().sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(default_params.population_priors_single().sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(default_params.germline_probability_mat_single().sum(), 16.0, 0.0001);
  BOOST_CHECK_CLOSE(default_params.germline_probability_mat().sum(), 256.0, 0.0001);
  // BOOST_CHECK_CLOSE(default_params.germline_probability_mat_num().sum(), 256.0, 0.0001);
  BOOST_CHECK_CLOSE(default_params.somatic_probability_mat().sum(), 16.0, 0.0001);
  // BOOST_CHECK_CLOSE(default_params.somatic_probability_mat_diag().sum(), 16.0, 0.0001);

  double theta = 0.05;
  double mutation_rate = 1e-6;
  double sequencing_error_rate = 0.25;
  RowVector4d frequencies_overloaded;
  frequencies_overloaded << 0.75, 0.0, 0.25, 0.0;
  MultinomialTrioModel overloaded_params(theta,
                                         mutation_rate,
                                         mutation_rate,
                                         sequencing_error_rate,
                                         frequencies_overloaded);

  double germline_mutation_rate_overloaded = overloaded_params.germline_mutation_rate();
  double exp_term_overloaded = exp(-4.0/3.0 * germline_mutation_rate_overloaded);
  double homozygous_match_overloaded = 0.25 + 0.75 * exp_term_overloaded;
  double heterozygous_match_overloaded = 0.25 + 0.25 * exp_term_overloaded;
  double mismatch_overloaded = 0.25 - 0.25 * exp_term_overloaded;

  BOOST_CHECK(overloaded_params.population_mutation_rate() == theta);
  BOOST_CHECK(germline_mutation_rate_overloaded == mutation_rate);
  BOOST_CHECK(overloaded_params.somatic_mutation_rate() == mutation_rate);
  BOOST_CHECK(overloaded_params.sequencing_error_rate() == sequencing_error_rate);
  BOOST_CHECK(overloaded_params.dirichlet_dispersion() == 1.0);
  BOOST_CHECK(overloaded_params.nucleotide_frequencies() == frequencies_overloaded);
  BOOST_CHECK(overloaded_params.homozygous_match() == homozygous_match_overloaded);
  BOOST_CHECK(overloaded_params.heterozygous_match() == heterozygous_match_overloaded);
  BOOST_CHECK(overloaded_params.mismatch() == mismatch_overloaded);
  BOOST_CHECK_CLOSE(overloaded_params.population_priors().sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(overloaded_params.population_priors_single().sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(overloaded_params.germline_probability_mat_single().sum(), 16.0, 0.0001);
  BOOST_CHECK_CLOSE(overloaded_params.germline_probability_mat().sum(), 256.0, 0.0001);
  // BOOST_CHECK_CLOSE(overloaded_params.germline_probability_mat_num().sum(), 256.0, 0.0001);
  BOOST_CHECK_CLOSE(overloaded_params.somatic_probability_mat().sum(), 16.0, 0.0001);
  // BOOST_CHECK_CLOSE(overloaded_params.somatic_probability_mat_diag().sum(), 16.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(TestMutationProbability) {
  MultinomialTrioModel params;
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  for (ReadDataVector data_vec : trio_vec) {
    double probability = params.MutationProbability(data_vec);
    BOOST_CHECK(probability >= 0.0 && probability <= 1.0);
  }
}

// Also tests read_dependent_data.cc
BOOST_AUTO_TEST_CASE(TestSetReadDependentData) {
  MultinomialTrioModel params;
  ReadDataVector data_vec = {{40, 0, 0, 0},
                             {40, 0, 0, 0},
                             {40, 0, 0, 0}};
  params.SetReadDependentData(data_vec);
  ReadDependentData data = params.read_dependent_data();
  BOOST_CHECK(EqualsReadDataVector(data_vec, data.read_data_vec));
  BOOST_CHECK(data.max_elements.size() > 0);
  BOOST_CHECK_CLOSE(data.sequencing_probability_mat.sum(), 3.0, 0.0001);
  BOOST_CHECK_CLOSE(data.child_somatic_probability.sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(data.mother_somatic_probability.sum(), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(data.father_somatic_probability.sum(), 1.0, 0.0001);
  BOOST_CHECK(data.denominator.child_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.denominator.mother_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.denominator.father_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.denominator.child_germline_probability.sum() > 0);
  BOOST_CHECK(data.denominator.parent_probability.sum() > 0);
  BOOST_CHECK(data.denominator.root_mat.sum() > 0);
  BOOST_CHECK(data.denominator.sum > 0.0);
  BOOST_CHECK(data.numerator.child_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.numerator.mother_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.numerator.father_zygotic_probability.sum() > 0);
  BOOST_CHECK(data.numerator.child_germline_probability.sum() > 0);
  BOOST_CHECK(data.numerator.parent_probability.sum() > 0);
  BOOST_CHECK(data.numerator.root_mat.sum() > 0);
  BOOST_CHECK(data.numerator.sum > 0.0);

  ReadDataVector data_vec2 = {{0, 0, 40, 0},
                              {40, 0, 0, 0},
                              {40, 0, 0, 0}};
  params.SetReadDependentData(data_vec2);
  BOOST_CHECK(!params.read_dependent_data().Equals(data));
}

BOOST_AUTO_TEST_CASE(TestGetAlphas) {
  MultinomialTrioModel params;
  Matrix16_4d test_alphas;
  Matrix16_4d alphas = params.alphas();
  double rate = params.sequencing_error_rate();
  double homozygous = 1.0 - rate;
  double mismatch = rate / 3.0;
  double heterozygous = 0.5 - mismatch;

  //             A             C             G             T
  test_alphas << homozygous,   mismatch,     mismatch,     mismatch,
                 heterozygous, heterozygous, mismatch,     mismatch,
                 heterozygous, mismatch,     heterozygous, mismatch,
                 heterozygous, mismatch,     mismatch,     heterozygous,

                 heterozygous, heterozygous, mismatch,     mismatch,
                 mismatch,     homozygous,   mismatch,     mismatch,
                 mismatch,     heterozygous, heterozygous, mismatch,
                 mismatch,     heterozygous, mismatch,     heterozygous,
            
                 heterozygous, mismatch,     heterozygous, mismatch,
                 mismatch,     heterozygous, heterozygous, mismatch,
                 mismatch,     mismatch,     homozygous,   mismatch,
                 mismatch,     mismatch,     heterozygous, heterozygous,

                 heterozygous, mismatch,     mismatch,     heterozygous,
                 mismatch,     heterozygous, mismatch,     heterozygous,
                 mismatch,     mismatch,     heterozygous, heterozygous,
                 mismatch,     mismatch,     mismatch,     homozygous;

  BOOST_CHECK(alphas == test_alphas);
}

BOOST_AUTO_TEST_CASE(TestEquals) {
  MultinomialTrioModel a;
  MultinomialTrioModel b(0.001, 2e-8, 2e-8, 0.005, {0.25, 0.25, 0.25, 0.25});
  MultinomialTrioModel c;
  c.set_population_mutation_rate(1e-6);
  BOOST_CHECK(a.Equals(b));
  BOOST_CHECK(!a.Equals(c));
}

BOOST_AUTO_TEST_CASE(test_set_population_mutation_rate) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  b.set_population_mutation_rate(1e-6);
  BOOST_CHECK(b.population_mutation_rate() == 1e-6);
  BOOST_CHECK(b.population_priors() != a.population_priors());
  BOOST_CHECK(b.population_priors_single() != a.population_priors_single());
}

BOOST_AUTO_TEST_CASE(test_set_germline_mutation_rate) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  b.set_germline_mutation_rate(1e-6);
  BOOST_CHECK(b.germline_mutation_rate() == 1e-6);
  BOOST_CHECK(b.homozygous_match() != a.homozygous_match());
  BOOST_CHECK(b.heterozygous_match() != a.heterozygous_match());
  BOOST_CHECK(b.mismatch() != a.mismatch());
  BOOST_CHECK(b.germline_probability_mat_single() != a.germline_probability_mat_single());
  BOOST_CHECK(b.germline_probability_mat() != a.germline_probability_mat());
  BOOST_CHECK(b.germline_probability_mat_num() != a.germline_probability_mat_num());
}

BOOST_AUTO_TEST_CASE(test_set_somatic_mutation_rate) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  b.set_somatic_mutation_rate(1e-6);
  BOOST_CHECK(b.somatic_mutation_rate() == 1e-6);
  BOOST_CHECK(b.somatic_probability_mat() != a.somatic_probability_mat());
  BOOST_CHECK(b.somatic_probability_mat_diag() != a.somatic_probability_mat_diag());
}

BOOST_AUTO_TEST_CASE(test_set_sequencing_error_rate) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  b.set_sequencing_error_rate(0.001);
  BOOST_CHECK(b.sequencing_error_rate() == 0.001);
  BOOST_CHECK(b.alphas() != a.alphas());
}

BOOST_AUTO_TEST_CASE(test_set_dirichlet_dispersion) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  b.set_dirichlet_dispersion(100000.0);
  BOOST_CHECK(b.dirichlet_dispersion() == 100000.0);
  // Uncomment this to check dirichlet multionomial on master branch.
  // BOOST_CHECK(b.alphas() != a.alphas());
}

BOOST_AUTO_TEST_CASE(test_set_nucleotide_frequencies) {
  MultinomialTrioModel a;
  MultinomialTrioModel b;
  RowVector4d frequencies;
  frequencies << 0.75, 0.0, 0.25, 0.0;
  b.set_nucleotide_frequencies(frequencies);
  BOOST_CHECK(b.nucleotide_frequencies() == frequencies);
  // Uncomment this to check dirichlet multionomial on master branch.
  // BOOST_CHECK(b.population_priors() != a.population_priors());
  // BOOST_CHECK(b.population_priors_single() != a.population_priors_single());
}
