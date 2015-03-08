/**
 * @file utility.h
 * @author Melissa Ip
 *
 * This file contains useful data structures, constants, and functions to
 * support the TrioModel including the Dirichlet multinomial.
 */
#ifndef UTILITY_H
#define UTILITY_H

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include <gsl/gsl_randist.h>  // Simulation and multinomial.

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

union ReadData {
  uint16_t reads[4];
  uint64_t key;
};

typedef Matrix<double, 1, 16> RowVector16d;
typedef Matrix<double, 1, 256> RowVector256d;
typedef Matrix<int, 16, 2, RowMajor> Matrix16_2i;
typedef Matrix<double, 3, 16, RowMajor> Matrix3_16d;
typedef Matrix<double, 4, 16, RowMajor> Matrix4_16d;
typedef Matrix<double, 16, 4, RowMajor> Matrix16_4d;
typedef Matrix<double, 16, 16, RowMajor> Matrix16_16d;
typedef Matrix<double, 16, 256, RowMajor> Matrix16_256d;
typedef Matrix<RowVector4d, 16, 16, RowMajor> Matrix16_16_4d;
typedef vector<ReadData> ReadDataVector;  // Contains child, mother, and father sequencing reads.
typedef vector<ReadDataVector> TrioVector;

// Forward declarations.
// Matrix16_2i GenotypeNumIndex();
uint16_t ToNucleotideIndex(char base);
Matrix16_16_4d ZeroMatrix16_16_4d();
void PrintMatrix16_16_4d(const Matrix16_16_4d &mat);
Matrix16_16_4d TwoParentCounts();
bool EqualsReadData(const ReadData &data1, const ReadData &data2);
bool EqualsReadDataVector(const ReadDataVector &data_vec1,
                          const ReadDataVector &data_vec2);
void PrintReadData(const ReadData &data);
void PrintReadDataVector(const ReadDataVector &data_vec);
MatrixXi EnumerateNucleotideCounts(int coverage);
ReadDataVector FourNucleotideCounts();
ReadDataVector GetPermutation(int arr[]);
ReadDataVector GetUniqueReadDataVector(const MatrixXi &mat);
int IndexOfReadDataVector(const ReadDataVector &data_vec, TrioVector trio_vec);
TrioVector GetTrioVector(int coverage);
bool HasZeroReadDataVector(const ReadDataVector &data_vec);
bool IsInVector(const RowVector4d &vec, double elem);
bool IsAlleleInParentGenotype(int child_nucleotide_idx, int parent_genotype_idx);
double DirichletMultinomialLog(const RowVector4d &alpha, const ReadData &data);
Matrix16_256d KroneckerProduct(const Matrix4_16d &mat);
Matrix16_16d KroneckerProduct(const Matrix4d &mat);
RowVector256d KroneckerProduct(const RowVector16d &vec1,
                               const RowVector16d &vec2);
bool Equal(double a, double b);
void Die(const char *msg);

/**
 * Global constants for specifying matrix size and iterating through numeric
 * representations of nucleotides and genotypes in lexicographical order.
 *
 * INDEX  GENOTYPE  NUCLEOTIDE 
 * 0      AA        A
 * 1      AC        C
 * 2      AG        G
 * 3      AT        T
 * 4      CA
 * 5      CC
 * 6      CG
 * 7      CT
 * 8      GA
 * 9      GC
 * 10     GG
 * 11     GT
 * 12     TA
 * 13     TC
 * 14     TG
 * 15     TT
 *
 * INDEX  GENOTYPE (paired parent genotype for use in 256 matrix)
 * 0      AAAA
 * 85     CCCC
 * 170    GGGG
 * 255    TTTT
 * i      Mother genotype*16 + Father genotype
 */
const int kAAAA = 0;
const int kCCCC = 85;
const int kGGGG = 170;
const int kTTTT = 255;
const int kNucleotideCount = 4;
const int kGenotypeCount = 16;
const int kGenotypePairCount = 256;
const int kTrioCount = 42875;
const double kThreshold = 0.01;  // Any greater probability than this number is printed.
const double kEpsilon = numeric_limits<double>::epsilon();
// const Matrix16_2i kGenotypeNumIndex = GenotypeNumIndex();
// const Matrix16_16_4d kTwoParentCounts = TwoParentCounts();

#endif
