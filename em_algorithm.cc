/**
 * @file em_algorithm.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm applied to a simplified and modified version of the trio model
 * (infinite sites model), which will be updated as necessary in Fall 2014 and
 * Spring 2015 to more complex and biologically realistic models including the
 * custom model using Dirichlet-multinomial approximations instead of
 * multinomial approximations.
 */
#include "trio_model.cc"  // FIXME: change to .h

// E-step methods.
double GetGermlineStatistic(TrioModel params);
Matrix16_256d GermlineMutationCounts(TrioModel params);
Matrix4_16d GermlineMutationCountsSingle(TrioModel params);
double GetSomaticStatistic(TrioModel params);
Matrix16_16d SomaticMutationCounts();


/**
 * Returns S_Germ the total number of nucleotide mismatches between parent and
 * child genotypes, or the expected number of germline mutations. This is the
 * sum of SM the number of nucleotide mismatches between m* and oa and SF the
 * number of nucleotide mismatches between f* and ob. Calculated during the
 * E-step of expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected germline mutations.
 */
double GetGermlineStatistic(TrioModel params) {
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_256d germline_mutation_counts = GermlineMutationCounts(params);
  Matrix16_256d germline_probability_mat = params.germline_probability_mat();
  
  RowVector256d s_germ = RowVector256d::Zero();

  double child_term1 = 0.0;
  double child_term2 = 0.0;

  // for each genotype x in population priors
  for (int x = 0; x < kGenotypePairCount; ++x) {
    // at top of branches
    // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = (germline_probability_mat(y, x) *  // germline genotype given parent pair=x
                     data->denominator.child_zygotic_probability(y));  // uses zygotic
      child_term2 = /* 0 + */ germline_mutation_counts(y, x);
      s_germ(x) += child_term1 * child_term2;
    }

    s_germ(x) /= data->denominator.child_germline_probability(x);
    
    // merges j branches
    // S(R_mom,R_dad,R_child, parent_pair=x)
    // mother and father have no expected germline mutations
    // s_germ(x) += 0.0;

    // at root of tree
    // S(R_mom,R_dad,R_child)
    s_germ(x) *= data->denominator.root_mat(x);  // root_mat includes population_priors
  }

  return s_germ.sum() / data->denominator.sum;
}

/**
 * Generates a 16 x 256 Eigen matrix holding the number of germline mutations
 * (0's, 0.5's, 1's, 2's). For example:
 *
 *           Parent pair
 * Germline  AAAA ACAA AGAA ATAA CAAA CCAA CGAA CTAA GAAA GCAA GGAA GTAA ...
 * AA        0
 * AC        
 * ...
 *
 * Assumes that the first genotype in the parent pair comes from the mother and
 * the first nucleotide of the child genotype comes from the mother.
 *
 * @param  params TrioModel object containing parameters.
 * @return        16 x 256 Eigen matrix holding the number of germline mutations.
 */
Matrix16_256d GermlineMutationCounts(TrioModel params) {
  Matrix4_16d germline_mutation_counts = GermlineMutationCountsSingle(params);
  Matrix16_256d mat = Matrix16_256d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypePairCount; ++j) {
      int child_allele1 = i / kNucleotideCount;
      int child_allele2 = i % kNucleotideCount;
      int mother_genotype = j / kGenotypeCount;
      int father_genotype = j % kGenotypeCount;
      mat(i, j) = (germline_mutation_counts(child_allele1, mother_genotype) +
                   germline_mutation_counts(child_allele2, father_genotype));
    }
  }

  return mat;
}

/**
 * Generates a 4 x 16 Eigen matrix holding the number of germline mutations
 * based on single parent.
 * 
 * Assumes that the first genotype in the parent pair comes from the mother and
 * the first nucleotide of the child genotype comes from the mother.
 *
 * @param  params TrioModel object containing parameters.
 * @return        4 x 16 Eigen matrix holding the number of germline mutations.
 */
Matrix4_16d GermlineMutationCountsSingle(TrioModel params) {
  Matrix4_16d mat = Matrix4_16d::Zero();
  double rate = params.no_match();
  Matrix4_16d germline_probability_mat = params.germline_probability_mat_single();

  for (int i = 0; i < kNucleotideCount; ++i) {  // child allele
    for (int j = 0; j < kGenotypeCount; ++j) {  // single parent genotype
      if (j % 5 == 0) {  // parent homozygous genotypes are divisible by 5
        if (IsAlleleInParentGenotype(i, j)) {
          mat(i, j) = 0.0;  // allele is inherited
        } else {
          mat(i, j) = 1.0;  // allele is mutated
        }
      } else {  // heterozygous parent genotypes
        mat(i, j) = rate * 0.5 / germline_probability_mat(i, j);
      }
    }
  }

  return mat;
}

/**
 * Returns S_Som the number of nucleotide mismatches between all x and x′,
 * or the expected number of somatic mutations. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected somatic mutations.
 */
double GetSomaticStatistic(TrioModel params) {
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_16d somatic_probability_mat = params.somatic_probability_mat();
  Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();
  Matrix16_256d germline_probability_mat = params.germline_probability_mat();

  RowVector256d s_som = RowVector256d::Zero();
  RowVector16d s_som_mother = RowVector16d::Zero();
  RowVector16d s_som_father = RowVector16d::Zero();
  RowVector16d s_som_child = RowVector16d::Zero();
  RowVector256d s_som_mother_x = RowVector256d::Zero();  // x refers to parent pair genotype
  RowVector256d s_som_father_x = RowVector256d::Zero();
  RowVector256d s_som_child_x = RowVector256d::Zero();

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // P(R|zygotic genotype)
  // S(R_mom, mom_zygotic=x), S(R_dad, dad_zygotic=x), S(R_child, child_zygotic=x)
  for (int x = 0; x < kGenotypeCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = somatic_probability_mat(x, y) * data->child_somatic_probability(y);
      mother_term1 = somatic_probability_mat(x, y) * data->mother_somatic_probability(y);
      father_term1 = somatic_probability_mat(x, y) * data->father_somatic_probability(y);

      child_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      mother_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      father_term2 = /* 0 + */ somatic_mutation_counts(x, y);

      s_som_child(x) += child_term1 * child_term2;  // sum over y_j
      s_som_mother(x) += mother_term1 * mother_term2;
      s_som_father(x) += father_term1 * father_term2;
    }

    s_som_child(x) /= data->denominator.child_zygotic_probability(x);
    s_som_mother(x) /= data->denominator.mother_zygotic_probability(x);
    s_som_father(x) /= data->denominator.father_zygotic_probability(x);
  }

  // for each genotype x in population priors
  for (int x = 0; x < kGenotypePairCount; ++x) {
    // at top of branches
    // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
    // where s_som_mother and s_som_father do not change
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = (germline_probability_mat(y, x) *  // germline genotype given parent pair=x
                     data->denominator.child_zygotic_probability(y));  // uses zygotic
      child_term2 = s_som_child(y) /* + 0 */;
      s_som_child_x(x) += child_term1 * child_term2;
    }

    s_som_child_x(x) /= data->denominator.child_germline_probability(x);
    
    // merges j branches
    // S(R_mom,R_dad,R_child, parent_pair=x)
    s_som(x) += (s_som_child_x(x) +
                 s_som_mother(x / kGenotypeCount) +
                 s_som_father(x % kGenotypeCount));

    // at root of tree
    // S(R_mom,R_dad,R_child)
    s_som(x) *= data->denominator.root_mat(x);  // root_mat includes population_priors
  }

  return s_som.sum() / data->denominator.sum;
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
Matrix16_16d SomaticMutationCounts() {
  Matrix16_16d mat = Matrix16_16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (i != j) {  // leaves diagonal as 0's
        int somatic_allele1 = i / kNucleotideCount;
        int somatic_allele2 = j % kNucleotideCount;
        int zygotic_allele1 = j / kNucleotideCount;
        int zygotic_allele2 = j % kNucleotideCount;

        // 1 mutation count per allele site segregation
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
