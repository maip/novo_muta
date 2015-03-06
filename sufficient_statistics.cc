/**
 * @file sufficient_statistics.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm.
 *
 * See top of sufficient_statistics.h for a complete description.
 */
#include "sufficient_statistics.h"


/**
 * Returns S_Theta the expected population mutation rate. Not listed in the
 * paper as a summary statistic. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param params TrioModel object containing parameters.
 * @return       Expected theta.
 */
double SufficientStatistics::GetPopulationMutationRateStatistic(const TrioModel &params) {
  const ReadDependentData data = params.read_dependent_data();
  const RowVector256d root_mat = data.denominator.root_mat;
  const double sum = data.denominator.sum;
  double AAAA = root_mat(0) / sum;
  double CCCC = root_mat(64) / sum;
  double GGGG = root_mat(128) / sum;
  double TTTT = root_mat(192) / sum;
  
  return AAAA + CCCC + GGGG + TTTT;
}

/**
 * Returns S_Het the number of nucleotide matches between a somatic heterozygous
 * genotype and its sequencing reads. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param params TrioModel object containing parameters.
 * @return       Number of expected heterozygous matches.
 */
double SufficientStatistics::GetHeterozygousStatistic(const TrioModel &params) {
  return GetSequencingErrorStatistic(params, 1);

}

/**
 * Returns S_Hom the number of nucleotide matches between a somatic homozygous
 * genotype and its sequencing reads. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param params TrioModel object containing parameters.
 * @return       Number of expected homozygous matches.
 */
double SufficientStatistics::GetHomozygousStatistic(const TrioModel &params) {
  return GetSequencingErrorStatistic(params, 0);

}

/**
 * Returns S_E the number of nucleotide mismatches between a somatic genotype
 * and its sequencing reads. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected mismatches.
 */
double SufficientStatistics::GetMismatchStatistic(const TrioModel &params) {
  return GetSequencingErrorStatistic(params, 2);
}

/**
 * Returns S_E, S_Hom, or S_Het based on parameters passed in.
 *
 * @param  params     TrioModel object containing parameters.
 * @para   match_flag Homozygous (0), or heterzygous (1) matches, or mismatches (2).
 * @return            Number of expected sequencing error statistic.
 */
double SufficientStatistics::GetSequencingErrorStatistic(const TrioModel &params,
                                                         int match_flag) {
  const ReadDependentData data = params.read_dependent_data();
  const Matrix16_16d somatic_probability_mat = params.somatic_probability_mat();
  const Matrix16_256d germline_probability_mat = params.germline_probability_mat();

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // S(R_mom, mom_somatic=x), S(R_dad, dad_somatic=x), S(R_child, child_somatic=x)
  RowVector256d s_e = RowVector256d::Zero();
  RowVector16d s_e_child = RowVector16d::Zero();
  RowVector16d s_e_mother = RowVector16d::Zero();
  RowVector16d s_e_father = RowVector16d::Zero();
  RowVector256d s_e_child_x = RowVector256d::Zero(); // Given x = parent pair genotype.
  
  Matrix3_16d matches;
  if (match_flag == 0) {
    matches = data.homozygous_matches;
  } else if (match_flag == 1) {
    matches = data.heterozygous_matches;
  } else if (match_flag == 2) {
    matches = data.mismatches;
  } else {
    cout << "Invalid matches." << endl;
    matches = Matrix3_16d::Zero();
  }

  // S(R_mom, mom_zygotic=x), S(R_dad, dad_zygotic=x), S(R_child, child_zygotic=x)
  for (int x = 0; x < kGenotypeCount; ++x) { // Zygotic genotype.
    for (int y = 0; y < kGenotypeCount; ++y) { // Somatic genotype.
      child_term1 = somatic_probability_mat(y, x) * data.child_somatic_probability(y);
      mother_term1 = somatic_probability_mat(y, x) * data.mother_somatic_probability(y);
      father_term1 = somatic_probability_mat(y, x) * data.father_somatic_probability(y);

      child_term2 = matches(0, y);
      mother_term2 = matches(1, y);
      father_term2 = matches(2, y);

      s_e_child(x) += child_term1 * child_term2; // Sums over y_j.
      s_e_mother(x) += mother_term1 * mother_term2;
      s_e_father(x) += father_term1 * father_term2;
    }

    // Dividing out child_zygotic_probability occurs at the next node.
    s_e_mother(x) /= data.denominator.mother_zygotic_probability(x);
    s_e_father(x) /= data.denominator.father_zygotic_probability(x);
  }

  // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
  // At top of branches where s_som_mother and s_som_father do not change,
  // because no somatic mutation or sequencing error can occur in germline.
  // For each parent pair genotype x in population priors.
  for (int x = 0; x < kGenotypePairCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) { // Child germline genotype.
      child_term1 = germline_probability_mat(y, x);
      child_term2 = s_e_child(y) /* + 0 */;
      s_e_child_x(x) += child_term1 * child_term2;
    }

    s_e_child_x(x) /= data.denominator.child_germline_probability(x); // Includes child_zygotic_probability.

    // S(R_mom,R_dad,R_child, parent_pair=x)
    // Merges j branches.
    s_e(x) += (s_e_child_x(x) +
               s_e_mother(x / kGenotypeCount) +
               s_e_father(x % kGenotypeCount));

    // S(R_mom,R_dad,R_child)
    // At root of tree.
    s_e(x) *= data.denominator.root_mat(x); // Includes population_priors.
  }

  return s_e.sum() / data.denominator.sum;
}

/**
    //Note: Many duplicate code here, GetMismatches(const ReadDataVector &data_vec) can call GetMismatches(const ReadData &data). Same for heter/homo
 * Returns S_Germ the total number of nucleotide mismatches between parent and
 * child genotypes, or the expected number of germline mutations. This is the
 * sum of SM the number of nucleotide mismatches between m* and oa and SF the
 * number of nucleotide mismatches between f* and ob. Calculated during the
 * E-step of expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected germline mutations.
 */
double SufficientStatistics::GetGermlineStatistic(const TrioModel &params) {
  const ReadDependentData data = params.read_dependent_data();
  const Matrix16_256d germline_mutation_counts = GermlineMutationCounts(params);
  const Matrix16_256d germline_probability_mat = params.germline_probability_mat();
  const RowVector256d population_priors = params.population_priors();
  RowVector256d s_germ = RowVector256d::Zero();

  double child_term1 = 0.0;
  double child_term2 = 0.0;
  
  // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
  // At top of branches.
  // For each parent pair genotype x in population priors.
  for (int x = 0; x < kGenotypePairCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) { // Child germline genotype.
      child_term1 = (germline_probability_mat(y, x) *
                     data.denominator.child_zygotic_probability(y));
      child_term2 = /* 0 + */ germline_mutation_counts(y, x);
      s_germ(x) += child_term1 * child_term2;
    }

    s_germ(x) *= data.denominator.mother_zygotic_probability(x / kGenotypeCount);
    s_germ(x) *= data.denominator.father_zygotic_probability(x % kGenotypeCount);
    s_germ(x) *= population_priors(x);
  }

  return s_germ.sum() / data.denominator.sum;
}

/**
 * Generates a 16 x 256 Eigen matrix holding the number of expected germline
 * mutations containing 0's, Het's, 1's, 1 + Het's, and 2's. The rows are child
 * germline genotypes and the columns are parent pair genotypes. Assumes that
 * the first genotype in the parent pair comes from the mother and the first
 * nucleotide of the child genotype comes from the mother, the second from the
 * father.
 *
 * @param  params TrioModel object containing parameters.
 * @return        16 x 256 Eigen matrix holding the number of germline mutations.
 */
Matrix16_256d SufficientStatistics::GermlineMutationCounts(const TrioModel &params) {
  Matrix16_256d mat = Matrix16_256d::Zero();
  Matrix4_16d germline_mutation_counts = GermlineMutationCountsSingle(params);
  int child_allele1 = 0;
  int child_allele2 = 0;
  int mother_genotype = 0;
  int father_genotype = 0;
  
  for (int x = 0; x < kGenotypeCount; ++x) {  // Child genotype.
    child_allele1 = x / kNucleotideCount;  // Inherited from mother.
    child_allele2 = x % kNucleotideCount;  // Inherited from father.

    for (int y = 0; y < kGenotypePairCount; ++y) {  // Parent pair genotype.
      mother_genotype = y / kGenotypeCount;
      father_genotype = y % kGenotypeCount;

      mat(x, y) = (germline_mutation_counts(child_allele1, mother_genotype) +
                   germline_mutation_counts(child_allele2, father_genotype));
    }
  }

  return mat;
}

/**
 * Generates a 4 x 16 Eigen matrix holding the number of expected germline
 * mutations or mismatches between the parent alleles and the zygote based on
 * a single parent. Assumes that the first genotype in the parent pair comes
 * from the mother and the first nucleotide of the child genotype comes from
 * the mother. For example:
 *
 *        Parent
 * Child  AA AC  AG  AT  CA  CC CG  CT  GA  GC  GG GT  TA  TC  TG  TT
 * A      0  Het Het Het Het 1  1   1   Het 1   1  1   Het 1   1   1
 * C      1  Het 1   1   Het 0  Het Het 1   Het 1  1   1   Het 1   1
 * G      1  1   Het 1   1   1  Het 1   Het Het 0  Het 1   1   Het 1
 * T      1  1   1   Het 1   1  1   Het 1   1   1  Het Het Het Het 0
 *
 * Het = 0.5 * no match / heterozygous match.
 *
 * @param  params TrioModel object containing parameters.
 * @return        4 x 16 Eigen matrix holding the number of germline mutations.
 */
Matrix4_16d SufficientStatistics::GermlineMutationCountsSingle(const TrioModel &params) {
  Matrix4_16d mat = Matrix4_16d::Zero();

  for (int i = 0; i < kNucleotideCount; ++i) {  // Child allele.
    for (int j = 0; j < kGenotypeCount; ++j) {  // Single parent genotype.
      if (j % 5 == 0) {  // Homozygous parent genotypes are divisible by 5.
        if (IsAlleleInParentGenotype(i, j)) {
          mat(i, j) = 0.0;  // Allele is inherited, ex: A in AA.
        } else {
          mat(i, j) = 1.0;  // Allele is mutated, ex: C in AA.
        }
      } else {  // Heterozygous parent genotype.
        if (IsAlleleInParentGenotype(i, j)) {
          // Allele is potentially mutated (Het), ex: C in AC.
          mat(i, j) = 0.5 * params.mismatch() / params.heterozygous_match();
        } else {
          mat(i, j) = 1.0;  // Allele is mutated, ex: A in CG.
        }
      }
    }
  }

  return mat;
}

/**
 * Returns S_Som, the number of nucleotide mismatches between all x and x′,
 * or the expected number of somatic mutations. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected somatic mutations.
 */
double SufficientStatistics::GetSomaticStatistic(const TrioModel &params) {
  const ReadDependentData data = params.read_dependent_data();
  const Matrix16_16d somatic_probability_mat = params.somatic_probability_mat();
  const Matrix16_256d germline_probability_mat = params.germline_probability_mat();
  const Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // S(R_mom, mom_somatic=x), S(R_dad, dad_somatic=x), S(R_child, child_somatic=x)
  RowVector256d s_e = RowVector256d::Zero();
  RowVector16d s_e_child = RowVector16d::Zero();
  RowVector16d s_e_mother = RowVector16d::Zero();
  RowVector16d s_e_father = RowVector16d::Zero();
  RowVector256d s_e_child_x = RowVector256d::Zero(); // Given x = parent pair genotype.
  
  // S(R_mom, mom_zygotic=x), S(R_dad, dad_zygotic=x), S(R_child, child_zygotic=x)
  for (int x = 0; x < kGenotypeCount; ++x) { // Zygotic genotype.
    for (int y = 0; y < kGenotypeCount; ++y) { // Somatic genotype.
      child_term1 = somatic_probability_mat(y, x) * data.child_somatic_probability(y);
      mother_term1 = somatic_probability_mat(y, x) * data.mother_somatic_probability(y);
      father_term1 = somatic_probability_mat(y, x) * data.father_somatic_probability(y);

      child_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      mother_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      father_term2 = /* 0 + */ somatic_mutation_counts(x, y);

      s_e_child(x) += child_term1 * child_term2; // Sums over y_j.
      s_e_mother(x) += mother_term1 * mother_term2;
      s_e_father(x) += father_term1 * father_term2;
    }

    // Dividing out child_zygotic_probability occurs at the next node.
    s_e_mother(x) /= data.denominator.mother_zygotic_probability(x);
    s_e_father(x) /= data.denominator.father_zygotic_probability(x);
  }

  // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
  // At top of branches where s_som_mother and s_som_father do not change,
  // because no somatic mutation or sequencing error can occur in germline.
  // For each parent pair genotype x in population priors.
  for (int x = 0; x < kGenotypePairCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) { // Child germline genotype.
      child_term1 = germline_probability_mat(y, x);
      child_term2 = s_e_child(y) /* + 0 */;
      s_e_child_x(x) += child_term1 * child_term2;
    }

    s_e_child_x(x) /= data.denominator.child_germline_probability(x); // Includes child_zygotic_probability.

    // S(R_mom,R_dad,R_child, parent_pair=x)
    // Merges j branches.
    s_e(x) += (s_e_child_x(x) +
               s_e_mother(x / kGenotypeCount) +
               s_e_father(x % kGenotypeCount));

    // S(R_mom,R_dad,R_child)
    // At root of tree.
    s_e(x) *= data.denominator.root_mat(x); // Includes population_priors.
  }

  return s_e.sum() / data.denominator.sum;
}

/**
 * Generates a 16 x 16 Eigen matrix holding the number of somatic mutations
 * where each allele site segregation is counted as a mutation in the infinite
 * sites model. For example:
 *
 *          Zygotic
 * Somatic  AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * AA       0  1  1  1  1  2  2  2  1  2  2  2  1  2  2  2
 * AC       1  0  1  1  2  1  2  2  2  1  2  2  2  1  2  2
 * ...            0
 *
 * @return  16 x 16 Eigen matrix holding the number of somatic mutations.
 */
Matrix16_16d SufficientStatistics::SomaticMutationCounts() {
  Matrix16_16d mat = Matrix16_16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (i != j) {  // Diagonal is filled with 0's for no mutation.
        int somatic_allele1 = i / kNucleotideCount;
        int somatic_allele2 = i % kNucleotideCount;
        int zygotic_allele1 = j / kNucleotideCount;
        int zygotic_allele2 = j % kNucleotideCount;

        // 1 mutation count per allele site segregation.
        if (somatic_allele1 != zygotic_allele1) {
          mat(i, j)++;
        }
        if (somatic_allele2 != zygotic_allele2) {
          mat(i, j)++;
        }
      }
    }
  }

  return mat;
}
