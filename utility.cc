/**
 * @file utility.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of utility.h.
 *
 * See top of utility.h for a complete description.
 */
#include "utility.h"


/**
 * Returns 16 x 2 Eigen matrix where the first dimension represents genotypes
 * and the second dimension represents the numeric index of each nucleotide in
 * that genotype.
 *
 * Used for generating the constant kGenotypeNumIndex only. This constant is
 * used to determine what nucleotides a genotype is composed of or what
 * genotype a pair of nucleotides is. The first allele can be determined by
 * integer division. The second allele can be determined by modulus.
 *
 * GENOTYPE  INDEX
 * AA        0, 0
 * AC        0, 1
 * AG        0, 2
 * AT        0, 3
 * ...       ...
 * TT        3, 3
 *
 * @return  16 x 2 Eigen matrix.
 */
// Matrix16_2i GenotypeNumIndex() {
//   Matrix16_2i mat;
//   //        A|    C|    G|    T|
//   mat << 0, 0, 0, 1, 0, 2, 0, 3,  // A
//          1, 0, 1, 1, 1, 2, 1, 3,  // C
//          2, 0, 2, 1, 2, 2, 2, 3,  // G
//          3, 0, 3, 1, 3, 2, 3, 3;  // T
//   return mat;
// }

/**
 * Returns numeric nucleotide index given a char. For use in VariantVisitor.
 *
 * @param  b Char.
 * @return   Numeric nucleotide index.
 */
uint16_t ToNucleotideIndex(char base) {
  switch (base) {        
  case 'A':
  case 'a':    
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '-':
  case 'N':
    return -1 ;
  default:  // Unknown base.
    return -1;
  }
}

/**
 * Returns a 16 x 16 x 4 Eigen matrix filled with zeros. The third dimension is
 * RowVector4d.
 *
 * @return  16 x 16 x 4 Eigen matrix filled with zeros.
 */
Matrix16_16_4d ZeroMatrix16_16_4d() {
  Matrix16_16_4d mat;
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      mat(i, j) = RowVector4d::Zero();
    }
  }
  return mat;
}

/**
 * Prints a Matrix16_16_4d object.
 *
 * @param  m Matrix16_16_4d object to be printed.
 */
void PrintMatrix16_16_4d(const Matrix16_16_4d &mat) {
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      cout << mat(i, j) << endl;
    }
  }
}

/**
 * Generates a 16 x 16 x 4 Eigen matrix. The first dimension represents the
 * mother genotype. The second dimension represents the father genotype. The
 * third dimension represents the nucleotide counts. For example, the (0, 0)
 * element represents AAAA, so its corresponding nucleotide counts would be
 * [4, 0, 0, 0].
 *
 * Used for generating the constant kTwoParentCounts only, which is used in
 * creating the population priors.
 *
 * @return  16 x 16 x 4 Eigen matrix of genotype counts where the (i, j) element
 *          is the count of nucleotide k of mother genotype i and father
 *          genotype j.
 */
Matrix16_16_4d TwoParentCounts() {
  Matrix16_16_4d counts = ZeroMatrix16_16_4d();
  RowVector4d vec = RowVector4d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      vec(i / kNucleotideCount)++;
      vec(i % kNucleotideCount)++;
      vec(j / kNucleotideCount)++;
      vec(j % kNucleotideCount)++;
      counts(i, j) = vec;
      vec = RowVector4d::Zero();
    }
  }
  return counts;
}

/**
 * Returns true if the two ReadDatas are equal. Assume the ReadData has a field
 * uint16_t reads[4].
 *
 * @param  data1 First ReadData.
 * @param  data2 Second ReadData.
 * @return       True if the ReadDatas are equal.
 */
bool EqualsReadData(const ReadData &data1, const ReadData &data2) {
  return equal(begin(data1.reads), end(data1.reads), begin(data2.reads));
}

/**
 * Returns true if the two ReadDataVectors are equal. Assume all ReadData has a
 * field uint16_t reads[4].
 *
 * @param  data_vec1 First ReadDataVector.
 * @param  data_vec2 Second ReadDataVector.
 * @return           True if the ReadDataVectors are equal.
 */
bool EqualsReadDataVector(const ReadDataVector &data_vec1,
                          const ReadDataVector &data_vec2) {
  bool table[3] = {EqualsReadData(data_vec1[0], data_vec2[0]),
                   EqualsReadData(data_vec1[1], data_vec2[1]),
                   EqualsReadData(data_vec1[2], data_vec2[2])};
  if (all_of(begin(table), end(table), [](bool i) { return i; })) {
    return true;
  } else {
    return false;
  }
}

/**
 * Prints ReadData with each read separated by spaces.
 *
 * @param  data ReadData to be printed.
 */
void PrintReadData(const ReadData &data) {
  cout << data.reads[0] << " " << data.reads[1] << " " << data.reads[2] << " "
       << data.reads[3] << endl;
}

/**
 * Prints ReadDataVector with each ReadData on a new line, and each read
 * separated by spaces.
 *
 * @param  data_vec ReadDataVector to be printed.
 */
void PrintReadDataVector(const ReadDataVector &data_vec) {
  for (const auto &data : data_vec) {
    PrintReadData(data);
  }
}

/**
 * Enumerates all possible nucleotide counts for an individual sequenced at
 * given coverage.
 *
 * The current implementation creates duplicate nucleotide counts.
 *
 * @param  coverage  Coverage or max nucleotide count. Must be at least 1.
 * @return           4^coverage x 4 Eigen matrix of nucleotide counts.
 */
MatrixXi EnumerateNucleotideCounts(int coverage) {
  Matrix4i identity_mat = Matrix4i::Identity();
  if (coverage == 1) {
    return identity_mat;
  } else {
    int rows = (int) pow(kNucleotideCount, coverage);
    MatrixXi counts(rows, kNucleotideCount);
    MatrixXi recursive = EnumerateNucleotideCounts(coverage - 1);
    for (int j = 0; j < recursive.rows(); ++j) {
      for (int i = 0; i < kNucleotideCount; ++i) {
        counts.row(i + j*kNucleotideCount) = (identity_mat.row(i) +
          recursive.row(j));
      }
    }
    return counts;
  }
}

/**
 * Returns all possible unique nucleotide counts at 4x coverage. The sum of each
 * nucleotide counts is 4.
 *
 * This function is an alternative to the current implementation of
 * EnumerateNucleotideCounts(), which has duplicate nucleotide counts and is
 * thus inefficient.
 *
 * @return  ReadDataVector containing unique ReadData counts at 4x coverage.
 */
ReadDataVector FourNucleotideCounts() {
  int a4_0[] = {4, 0, 0, 0};
  int a3_1[] = {3, 1, 0, 0};
  int a2_2[] = {2, 2, 0, 0};
  int a2_1_1[] = {2, 1, 1, 0};

  ReadDataVector a4_0_counts = GetPermutation(a4_0);
  ReadDataVector a3_1_counts = GetPermutation(a3_1);
  ReadDataVector a2_2_counts = GetPermutation(a2_2);
  ReadDataVector a2_1_1_counts = GetPermutation(a2_1_1);

  ReadDataVector vec;
  vec.reserve(a4_0_counts.size() + a3_1_counts.size() +
              a2_2_counts.size() + a2_1_1_counts.size() + 1);
  vec.insert(vec.end(), a4_0_counts.begin(), a4_0_counts.end());
  vec.insert(vec.end(), a3_1_counts.begin(), a3_1_counts.end());
  vec.insert(vec.end(), a2_2_counts.begin(), a2_2_counts.end());
  vec.insert(vec.end(), a2_1_1_counts.begin(), a2_1_1_counts.end());

  ReadData a1_1_1_1 = {1, 1, 1, 1};  // Does not need to be sorted.
  vec.push_back(a1_1_1_1);
  return vec;
}

/**
 * Creates a ReadDataVector holding all unique combinations by rearranging a
 * given array of 4 integers.
 *
 * @param  arr Integer array holding nucleotide counts.
 * @return     ReadDataVector with all unique combinations of arr.
 */
ReadDataVector GetPermutation(int arr[]) {
  ReadData data;
  ReadDataVector vec;
  sort(arr, arr + kNucleotideCount);

  do {
    data.reads[0] = arr[0];
    data.reads[1] = arr[1];
    data.reads[2] = arr[2];
    data.reads[3] = arr[3];
    vec.push_back(data);
  } while (next_permutation(arr, arr + kNucleotideCount));

  return vec;
}

/**
 * Converts nucleotide counts in mat to ReadData and appends to ReadDataVector
 * if it does not already exists in the ReadDataVector.
 *
 * @param  mat Eigen matrix containing nucleotide counts.
 * @return     ReadDataVector containing all unique nucleotide counts converted
 *             to ReadData.
 */
ReadDataVector GetUniqueReadDataVector(const MatrixXi &mat) {
  ReadDataVector data_vec;
  for (int i = 0; i < mat.rows(); ++i) {
    bool is_in_vec = false;
    for (auto data : data_vec) {
      if (data.reads[0] == mat(i, 0) && data.reads[1] == mat(i, 1) && 
          data.reads[2] == mat(i, 2) && data.reads[3] == mat(i, 3)) {
        is_in_vec = true;
      }
    }

    if (!is_in_vec) {
      ReadData data = {0};
      for (int j = 0; j < kNucleotideCount; ++j) {
        data.reads[j] = mat(i, j);
      }
      data_vec.push_back(data);
    }
  }
  return data_vec;
}

/**
 * Returns the index of a ReadDataVector in the TrioVector of all trios at 4x
 * coverage. Assumes the ReadDataVector has 4x coverage.
 *
 * @param  data_vec ReadDataVector.
 * @param  trio_vec TrioVector to look through.
 * @return          Index of ReadDataVector in TrioVector list.    
 */
int IndexOfReadDataVector(const ReadDataVector &data_vec, TrioVector trio_vec) {
  for (int i = 0; i < trio_vec.size(); ++i) {
    if (EqualsReadDataVector(data_vec, trio_vec[i])) {
      return i;
    }
  }
  return -1;
}

/**
 * Returns all possible and unique trio sets of sequencing counts for an
 * individual sequenced at given coverage. Assume  coverage is at least 1,
 * usually 4.
 *
 * @param  coverage Coverage or max nucleotide count.
 * @return          Vector of ReadDataVector.
 */
TrioVector GetTrioVector(int coverage) {
  TrioVector trio_vec;
  ReadDataVector data_vec;

  if (coverage == kNucleotideCount) {
    data_vec = FourNucleotideCounts();
  } else {  // Less efficient.
    MatrixXi mat = EnumerateNucleotideCounts(coverage);
    data_vec = GetUniqueReadDataVector(mat);
  }
  
  for (ReadData data1 : data_vec) {
    for (ReadData data2 : data_vec) {
      for (ReadData data3 : data_vec) {
        trio_vec.push_back({ data1, data2, data3 });
      }
    }
  }
  return trio_vec;
}

/**
 * Returns true if ReadDataVector contains a zero vector.
 *
 * @param  data_vec ReadDataVector
 * @return          True if ReadDataVector contains a zero vector.
 */
bool HasZeroReadDataVector(const ReadDataVector &data_vec) {
  ReadData zero_data = {0, 0, 0, 0};
  for (auto data : data_vec) {
    if (EqualsReadData(data, zero_data)) {
      return true;
    }
  }
  return false;
}

/**
 * Returns true if the element is in the RowVector.
 *
 * @param  vec  Eigen RowVector. 
 * @param  elem Element to look for in RowVector.
 * @return      True if element is in the RowVector.
 */
bool IsInVector(const RowVector4d &vec, double elem) {
  for (int i = 0; i < kNucleotideCount; ++i) {
    if (vec(i) == elem) {  // elem is [2.0, 4.0]
      return true;
    }
  }
  return false;
}

/**
 * Checks if the child nucleotide is in the parent genotype.
 *
 * @param  child_nucleotide_idx Index of child allele.
 * @param  parent_genotype_idx  Index of parent genotype.
 * @return                      True if child allele is in parent genotype.
 */
bool IsAlleleInParentGenotype(int child_nucleotide_idx, int parent_genotype_idx) {
  int parent_allele1 = parent_genotype_idx / kNucleotideCount;
  int parent_allele2 = parent_genotype_idx % kNucleotideCount;
  if (child_nucleotide_idx == parent_allele1 ||
      child_nucleotide_idx == parent_allele2) {
    return true;
  } else {
    return false;
  }
}

/**
 * Calculates probability from the probability density function (pdf):
 *
 * \frac{\gamma(\theta)}{\gamma(\theta + N)} *
 *     \Pi_{i = A, C, G, T} \frac{\gamma(\alpha_i * \theta + n_i)}
 *                               {\gamma(\alpha_i * \theta}
 *
 * @param  alpha RowVector of doubles representing frequencies for each
 *               category in the Dirichlet multinomial.
 * @param  data  Read counts or samples for each category in the multinomial.
 * @return       log_e(P) where P is the value calculated from the pdf.
 */
double DirichletMultinomialLog(const RowVector4d &alpha, const ReadData &data) {
  double a = alpha.sum();
  int n = data.reads[0] + data.reads[1] + data.reads[2] + data.reads[3];
  double constant_term = lgamma(a) - lgamma(n + a);
  double product_term = 0.0;
  for (int i = 0; i < kNucleotideCount; ++i) {
    product_term += lgamma(alpha(i) + data.reads[i]) - lgamma(alpha(i));
  }
  return constant_term + product_term;
}

/**
 * Calculates the Kronecker product of a matrix with itself.
 *
 * Matrix sizes are specific to germline_probability_mat_.
 *
 * @param  arr 4 x 16 Eigen matrix.
 * @return     16 x 256 Eigen matrix.
 */
Matrix16_256d KroneckerProduct(const Matrix4_16d &mat) {
  Matrix16_256d kronecker_product;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 16; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 16; ++l) {
          kronecker_product(i*4 + k, j*16 + l) = mat(i, j) * mat(k, l);
        }
      }
    }
  }
  return kronecker_product;
}

/**
 * Calculates the Kronecker product of a matrix with itself.
 *
 * Matrix sizes are specific to somatic_probability_mat_.
 *
 * @param  arr 4 x 4 Eigen matrix.
 * @return     16 x 16 Eigen matrix.
 */
Matrix16_16d KroneckerProduct(const Matrix4d &mat) {
  Matrix16_16d kronecker_product;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 4; ++l) {
          kronecker_product(i*4 + k, j*4 + l) = mat(i, j) * mat(k, l);
        }
      }
    }
  }
  return kronecker_product;
}

/**
 * Calculates the Kronecker product of a RowVector with itself.
 *
 * Matrix sizes are specific to the parent probability matrix.
 *
 * @param  arr1 1 x 16 Eigen RowVector.
 * @param  arr2 1 x 16 Eigen RowVector.
 * @return      1 x 256 Eigen RowVector.
 */
RowVector256d KroneckerProduct(const RowVector16d &vec1,
                               const RowVector16d &vec2) {
  RowVector256d kronecker_product;
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      kronecker_product(i*16 + j) = vec1(i) * vec2(j);
    }
  }
  return kronecker_product;
}

/**
 * Returns true if the two given doubles are equal to each other within epsilon
 * precision.
 *
 * @param  a First double.
 * @param  b Second double.
 * @return   True if the two doubles are equal to each other.
 */
bool Equal(double a, double b) {
  return fabs(a - b) < kEpsilon;
}

/**
 * Prints user defined errors and exits the program.
 *
 * @param  msg Message to be printed.
 */
void Die(const char *msg) {
  if (errno == EDOM) {
    perror(msg);
  } else {
    printf("ERROR: %s\n", msg);
  }
  exit(EXIT_FAILURE);
}
