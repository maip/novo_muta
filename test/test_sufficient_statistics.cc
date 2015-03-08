/**
 * @file test_sufficient_statistics.cc
 * @author Melissa Ip
 *
 * This file tests the functions in sufficient_statistics.h.
 */
#define BOOST_TEST_MODULE TestSufficientStatistics

#include <boost/test/unit_test.hpp>

#include "sufficient_statistics.h"


BOOST_AUTO_TEST_CASE(TestGetThetaStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_theta = SufficientStatistics::GetThetaStatistic(params);
  BOOST_CHECK(s_theta >= 0.0);
}

BOOST_AUTO_TEST_CASE(TestGetHeterozygousStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_het = SufficientStatistics::GetHeterozygousStatistic(params);
  BOOST_CHECK(s_het >= 0.0);
}

BOOST_AUTO_TEST_CASE(TestGetHomozygousStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_hom = SufficientStatistics::GetHomozygousStatistic(params);
  BOOST_CHECK(s_hom >= 0.0);
}

BOOST_AUTO_TEST_CASE(TestGetMismatchStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_e = SufficientStatistics::GetMismatchStatistic(params);
  BOOST_CHECK(s_e >= 0.0);
}

// BOOST_AUTO_TEST_CASE(TestGetSequencingErrorStatistic) {}

BOOST_AUTO_TEST_CASE(TestGetGermlineStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_germ = SufficientStatistics::GetGermlineStatistic(params);
  BOOST_CHECK(s_germ >= 0.0);
}

BOOST_AUTO_TEST_CASE(TestGermlineMutationCounts) {
  TrioModel params;
  Matrix16_256d counts = SufficientStatistics::GermlineMutationCounts(params);
  Matrix16_256d test_counts = Matrix16_256d::Zero();
  int child_allele1 = 0;
  int child_allele2 = 0;
  int mother_genotype = 0;
  int father_genotype = 0;

  for (int x = 0; x < kGenotypeCount; ++x) {
    child_allele1 = x / kNucleotideCount;
    child_allele2 = x % kNucleotideCount;

    for (int y = 0; y < kGenotypePairCount; ++y) {
      mother_genotype = y / kGenotypeCount;
      father_genotype = y % kGenotypeCount;

      test_counts(x, y) = (counts(child_allele1, mother_genotype) +
                           counts(child_allele2, father_genotype));
    }
  }
  BOOST_CHECK(counts == test_counts);
}

BOOST_AUTO_TEST_CASE(TestGermlineMutationCountsSingle) {
  TrioModel params;
  double heterozygous = 0.5 * params.mismatch() / params.heterozygous_match();
  Matrix4_16d test_counts = Matrix4_16d::Zero();
  test_counts << 0.0,          heterozygous, heterozygous, heterozygous,  // A
                 heterozygous, 1.0,          1.0,          1.0,
                 heterozygous, 1.0,          1.0,          1.0,
                 heterozygous, 1.0,          1.0,          1.0,
                 1.0,          heterozygous, 1.0,          1.0,           // C
                 heterozygous, 0.0,          heterozygous, heterozygous,
                 1.0,          heterozygous, 1.0,          1.0,
                 1.0,          heterozygous, 1.0,          1.0,
                 1.0,          1.0,          heterozygous, 1.0,           // G
                 1.0,          1.0,          heterozygous, 1.0,
                 heterozygous, heterozygous, 0.0,          heterozygous,
                 1.0,          1.0,          heterozygous, 1.0,
                 1.0,          1.0,          1.0,          heterozygous,  // T
                 1.0,          1.0,          1.0,          heterozygous,
                 1.0,          1.0,          1.0,          heterozygous,
                 heterozygous, heterozygous, heterozygous, 0.0;
  Matrix4_16d counts = SufficientStatistics::GermlineMutationCountsSingle(params);
  BOOST_CHECK(counts == test_counts);
}

BOOST_AUTO_TEST_CASE(TestGetSomaticStatistic) {
  TrioModel params;
  ReadDataVector data = {{20, 0, 20, 0},
                         {40, 0, 0,  0},
                         {40, 0, 0,  0}};
  params.SetReadDependentData(data);
  double s_som = SufficientStatistics::GetSomaticStatistic(params);
  BOOST_CHECK(s_som >= 0.0);
}

BOOST_AUTO_TEST_CASE(TestSomaticMutationCounts) {
  Matrix16_16d counts = SufficientStatistics::SomaticMutationCounts();
  BOOST_CHECK(counts.diagonal().sum() == 0.0);
  int count = 0;

  for (int i = 0; i < kGenotypeCount; ++i) {
    BOOST_CHECK(counts.row(i).sum() == 24.0);
    BOOST_CHECK(counts.col(i).sum() == 24.0);
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (i != j) {
        int somatic_allele1 = i / kNucleotideCount;
        int somatic_allele2 = i % kNucleotideCount;
        int zygotic_allele1 = j / kNucleotideCount;
        int zygotic_allele2 = j % kNucleotideCount;

        if (somatic_allele1 != zygotic_allele1) {
          count++;
        }
        if (somatic_allele2 != zygotic_allele2) {
          count++;
        }

        BOOST_CHECK(counts(i, j) == count);
        count = 0;  // Resets count for next position.
      }
    }
  }
}