/**
 * @file test_utility.cc
 * @author Melissa Ip
 *
 * This file tests the functions in utility.h.
 */
#define BOOST_TEST_MODULE TestUtility

#include <boost/test/unit_test.hpp>

#include "utility.h"


// BOOST_AUTO_TEST_CASE(TestGenotypeNumIndex) {
//   Matrix16_2i mat = GenotypeNumIndex();
//   for (int i = 0; i < kNucleotideCount; ++i) {
//     BOOST_CHECK(mat(i, 0) == 0);
//     BOOST_CHECK(mat(i+4, 0) == 1);
//     BOOST_CHECK(mat(i+8, 0) == 2);
//     BOOST_CHECK(mat(i+12, 0) == 3);
//     BOOST_CHECK(mat(i, 1) == i);
//     BOOST_CHECK(mat(i+4, 1) == i);
//     BOOST_CHECK(mat(i+8, 1) == i);
//     BOOST_CHECK(mat(i+12, 1) == i);
//   }
// }

BOOST_AUTO_TEST_CASE(TestZeroMatrix16_16_4d) {
  Matrix16_16_4d mat = ZeroMatrix16_16_4d();
  RowVector4d vec = RowVector4d::Zero();
  BOOST_REQUIRE(mat.rows() == kGenotypeCount);
  BOOST_REQUIRE(mat.cols() == kGenotypeCount);
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      RowVector4d nested = mat(i, j);
      BOOST_CHECK(nested.size() == kNucleotideCount);
      BOOST_CHECK(nested == vec);
    }
  }
}

// BOOST_AUTO_TEST_CASE(TestPrintMatrix16_16_4d) {}

BOOST_AUTO_TEST_CASE(TestTwoParentCounts) {
  Matrix16_16_4d mat = TwoParentCounts();
  RowVector4d vec = RowVector4d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      vec(i / kNucleotideCount)++;
      vec(i % kNucleotideCount)++;
      vec(j / kNucleotideCount)++;
      vec(j % kNucleotideCount)++;
      RowVector4d nested = mat(i, j);
      BOOST_CHECK(nested == vec);
      vec = RowVector4d::Zero();
    }
  }

  RowVector4d expected_4_ACGT[] = { RowVector4d {4,0,0,0},
                                    RowVector4d {0,4,0,0},
                                    RowVector4d {0,0,4,0},
                                    RowVector4d {0,0,0,4}};
  BOOST_CHECK(mat(0) == expected_4_ACGT[0]);
  BOOST_CHECK(mat(85) == expected_4_ACGT[1]);
  BOOST_CHECK(mat(170) == expected_4_ACGT[2]);
  BOOST_CHECK(mat(255) == expected_4_ACGT[3]);
}

BOOST_AUTO_TEST_CASE(TestEqualsReadData) {
  ReadData a = {1, 1, 1, 1};
  ReadData b = {1, 1, 1, 1};
  ReadData c = {2, 2, 2, 2};
  BOOST_CHECK(EqualsReadData(a, b));
  BOOST_CHECK(!EqualsReadData(a, c));
}

BOOST_AUTO_TEST_CASE(TestEqualsReadDataVector) {
  ReadDataVector a = {{1, 1, 1, 1},
                      {1, 1, 1, 1},
                      {1, 1, 1, 1}};
  ReadDataVector b = {{1, 1, 1, 1},
                      {1, 1, 1, 1},
                      {1, 1, 1, 1}};
  ReadDataVector c = {{2, 2, 2, 2},
                      {2, 2, 2, 2},
                      {2, 2, 2, 2}};
  BOOST_CHECK(EqualsReadDataVector(a, b));
  BOOST_CHECK(!EqualsReadDataVector(a, c));
}

// BOOST_AUTO_TEST_CASE(TestPrintReadData) {}

// BOOST_AUTO_TEST_CASE(TestPrintReadDataVector) {}

BOOST_AUTO_TEST_CASE(TestEnumerateNucleotideCounts) {
  MatrixXi mat = EnumerateNucleotideCounts(kNucleotideCount);
  for (int i = 0; i < mat.rows(); ++i) {
    BOOST_CHECK(mat.row(i).sum() == kNucleotideCount);
  }
}

BOOST_AUTO_TEST_CASE(TestFourNucleotideCounts) {
  ReadDataVector vec = FourNucleotideCounts();
  BOOST_REQUIRE(vec.size() == 35);
  for (ReadData data : vec) {
    BOOST_CHECK(data.reads[0] + data.reads[1] +
                data.reads[2] + data.reads[3] == kNucleotideCount);
  }
}

BOOST_AUTO_TEST_CASE(TestGetTrioVector) {
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  BOOST_REQUIRE(trio_vec.size() == kTrioCount);
  ReadDataVector vec = {{0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0}};
  BOOST_CHECK(IndexOfReadDataVector(vec, trio_vec) == -1);
}

BOOST_AUTO_TEST_CASE(TestGetUniqueReadDataVector) {
  MatrixXi mat(4, 4);
  mat << 0, 0, 0, 0,
         1, 1, 1, 1,
         0, 0, 0, 0,
         0, 1, 0, 1;
  ReadDataVector vec = GetUniqueReadDataVector(mat);
  BOOST_CHECK(vec.size() == 3);
}

BOOST_AUTO_TEST_CASE(TestIndexOfReadDataVector) {
  ReadDataVector a = {{0, 0, 0, 0},
                      {0, 0, 0, 0},
                      {0, 0, 0, 0}};
  ReadDataVector b = {{1, 1, 1, 1},
                      {1, 1, 1, 1},
                      {1, 1, 1, 1}};
  ReadDataVector c = {{2, 2, 2, 2},
                      {2, 2, 2, 2},
                      {2, 2, 2, 2}};
  TrioVector trio_vec = {a, b};
  BOOST_CHECK(IndexOfReadDataVector(b, trio_vec) == 1);
  BOOST_CHECK(IndexOfReadDataVector(c, trio_vec) == -1);
}

BOOST_AUTO_TEST_CASE(TestIsInVector) {
  RowVector4d vec;
  vec << 0.0, 1.0, 2.0, 3.0;
  BOOST_CHECK(IsInVector(vec, 1.0));
  BOOST_CHECK(!IsInVector(vec, 4.0));
}

BOOST_AUTO_TEST_CASE(TestIsAlleleInParentGenotype) {
  for (int i = 0; i < kGenotypeCount; ++i) {
    BOOST_CHECK(IsAlleleInParentGenotype(i / kNucleotideCount, i));
    BOOST_CHECK(IsAlleleInParentGenotype(i % kNucleotideCount, i));
  }
}

BOOST_AUTO_TEST_CASE(TestDirichletMultinomialLog) {
  RowVector4d alpha;
  alpha << 0.25, 0.25, 0.25, 0.25;
  ReadDataVector vec = FourNucleotideCounts();
  double probability = 0.0;
  for (ReadData data : vec) {
    probability = exp(DirichletMultinomialLog(alpha, data));
    BOOST_CHECK(probability >= 0.0 && probability <= 1.0);
  }
}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct) {
  Matrix4_16d mat;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 16; ++j) {
      mat(i, j) = j;
    }
  }
  Matrix16_256d kron = KroneckerProduct(mat);
  BOOST_REQUIRE(kron.rows() == 16);
  BOOST_REQUIRE(kron.cols() == 256);
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 256; ++j) {
      BOOST_CHECK(kron(i, j) == mat(i/4, j/16) * mat(i%4, j%16));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct2) {
  Matrix4d mat;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat(i, j) = j;
    }
  }
  Matrix16_16d kron = KroneckerProduct(mat);
  BOOST_REQUIRE(kron.rows() == 16);
  BOOST_REQUIRE(kron.cols() == 16);
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      BOOST_CHECK(kron(i, j) == mat(i/4, j/4) * mat(i%4, j%4));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct3) {
  RowVector16d vec;
  for (int i = 0; i < 16; ++i) {
    vec(i) = i;
  }
  RowVector256d kron = KroneckerProduct(vec, vec);
  BOOST_REQUIRE(kron.rows() == 1);
  BOOST_REQUIRE(kron.cols() == 256);
  for (int i = 0; i < 256; ++i) {
    BOOST_CHECK(kron(i) == vec(i/16) * vec(i%16));
  }
}

BOOST_AUTO_TEST_CASE(TestEqual) {
  BOOST_CHECK(Equal(0.0, 0.0));
  BOOST_CHECK(!Equal(0.0, 1.0));
  BOOST_CHECK(Equal(0.0, 2.22e-16));
  BOOST_CHECK(!Equal(0.0, 2.23e-16));
}

// BOOST_AUTO_TEST_CASE(TestDie) {}