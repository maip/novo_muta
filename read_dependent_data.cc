/**
 * @file read_dependent_data.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the ReadDependentData class.
 *
 * See top of read_dependent_data.h for a complete description.
 */
#include "read_dependent_data.h"


/**
 * Default constructor.
 */
ReadDependentData::ReadDependentData() {
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
  homozygous_matches = Matrix3_16d::Zero();
  heterozygous_matches = Matrix3_16d::Zero();
  mismatches = Matrix3_16d::Zero();
}

/**
 * Constructor takes in ReadDataVector.
 */
ReadDependentData::ReadDependentData(const ReadDataVector &data_vec)
    : read_data_vec{data_vec} {
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
  homozygous_matches = GetHomozygousMatches();
  heterozygous_matches = GetHeterozygousMatches();
  mismatches = GetMismatches();
};

/**
 * Sums all nucleotide counts in each ReadData and subtracts out the number of
 * nucleotides that match the genotype. It does not subtract twice for
 * homozygous genotypes. Returns 3 x 16 Eigen matrix holding number of
 * mismatches per genotype for each read. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 11 1  11 10 1  21 21 20 11 21 31 30 10 20 30 30
 *
 * @return  3 x 16 Eigen matrix containing number of mismatches per genotype
 *          for each read data.
 */
Matrix3_16d ReadDependentData::GetMismatches() {
  Matrix3_16d s_e = Matrix3_16d::Zero();
  for (int i = 0; i < 3; ++i) {
    ReadData data = read_data_vec[i];
    for (int j = 0; j < kGenotypeCount; ++j) {
      int allele1 = j / kNucleotideCount;
      int allele2 = j % kNucleotideCount;

      s_e(i, j) += data.reads[0] + data.reads[1] + data.reads[2] + data.reads[3];
      s_e(i, j) -= data.reads[allele1];  // Homozygous.
      if (allele1 != allele2) {
        s_e(i, j) -= data.reads[allele2];  // Hetereogyzous.
      }
    }
  }

  return s_e;
}

/**
 * Returns 3 x 16 Eigen matrix holding number of heterozygous matches per
 * genotype for each read. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 0  30 20 21 30 0  10 11 20 10 0  1  21 11 1  0 
 *
 * @return  3 x 16 Eigen matrix containing number of heterozygous matches per
 *          genotype for each read data.
 */
Matrix3_16d ReadDependentData::GetHeterozygousMatches() {
  Matrix3_16d s_het = Matrix3_16d::Zero();
  for (int i = 0; i < 3; ++i) {
    ReadData data = read_data_vec[i];
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (j % 5 != 0) {  // Heterozygous genotypes are not divisible by 5.
        int allele1 = j / kNucleotideCount;
        int allele2 = j % kNucleotideCount;
        s_het(i, j) += data.reads[allele1] + data.reads[allele2];
      }
    }
  }

  return s_het;
}

/**
 * Returns 3 x 16 Eigen matrix holding number of homozygous matches per
 * genotype for each read. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 20 0  0  0  0  10 0  0  0  0  0  0  0  0  0  1 
 *
 * @return  3 x 16 Eigen matrix containing number of homozygous matches per
 *          genotype for each read data.
 */
Matrix3_16d ReadDependentData::GetHomozygousMatches() {
  Matrix3_16d s_hom = Matrix3_16d::Zero();
  for (int i = 0; i < 3; ++i) {
    ReadData data = read_data_vec[i];
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (j % 5 == 0) {  // Homozygous genotypes are divisible by 5.
        s_hom(i, j) += data.reads[j / kNucleotideCount];
      }
    }
  }

  return s_hom;
}

/**
 * Returns true if the two ReadDependentData objects are equal to each other.
 *
 * @param other ReadDependentData object to be compared.
 * @return True if the two ReadDependentData objects are equal to each other.
 */
bool ReadDependentData::Equals(const ReadDependentData &other) {
  bool attr_table[20] = {
    EqualsReadDataVector(read_data_vec, other.read_data_vec),
    max_elements == other.max_elements,  // Better comparison of vectors?
    sequencing_probability_mat.isApprox(other.sequencing_probability_mat, kEpsilon),
    child_somatic_probability.isApprox(other.child_somatic_probability, kEpsilon),
    mother_somatic_probability.isApprox(other.mother_somatic_probability, kEpsilon),
    father_somatic_probability.isApprox(other.father_somatic_probability, kEpsilon),
    denominator.child_zygotic_probability.isApprox(other.denominator.child_zygotic_probability, kEpsilon),
    denominator.mother_zygotic_probability.isApprox(other.denominator.mother_zygotic_probability, kEpsilon),
    denominator.father_zygotic_probability.isApprox(other.denominator.father_zygotic_probability, kEpsilon),
    denominator.child_germline_probability.isApprox(other.denominator.child_germline_probability, kEpsilon),
    denominator.parent_probability.isApprox(other.denominator.parent_probability, kEpsilon),
    denominator.root_mat.isApprox(other.denominator.root_mat, kEpsilon),
    Equal(denominator.sum, other.denominator.sum),
    numerator.child_zygotic_probability.isApprox(other.numerator.child_zygotic_probability, kEpsilon),
    numerator.mother_zygotic_probability.isApprox(other.numerator.mother_zygotic_probability, kEpsilon),
    numerator.father_zygotic_probability.isApprox(other.numerator.father_zygotic_probability, kEpsilon),
    numerator.child_germline_probability.isApprox(other.numerator.child_germline_probability, kEpsilon),
    numerator.parent_probability.isApprox(other.numerator.parent_probability, kEpsilon),
    numerator.root_mat.isApprox(other.numerator.root_mat, kEpsilon),
    Equal(numerator.sum, other.numerator.sum)
  };

  if (all_of(begin(attr_table), end(attr_table), [](bool i) { return i; })) {
    return true;
  } else {
    return false;
  }
}