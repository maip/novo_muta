/**
 * @file multinomial_trio_model.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the MultinomialTrioModel class.
 *
 * See top of multinomial_trio_model.h for a complete description.
 */
#include "multinomial_trio_model.h"


/**
 * Default constructor.
 * 
 * sequencing_probability_mat is created or updated if sequencing_error_rate_
 * or dirichlet_dispersion_ is changed when MutationProbability() or
 * SetReadDependentData() is called. dirichlet_dispersion_ is not used in the
 * infinite sites model version.
 */
MultinomialTrioModel::MultinomialTrioModel()
    : TrioModel() {
  set_dirichlet_dispersion(1.0);
  population_priors_ = TrioModel::PopulationPriors();
  population_priors_single_ = TrioModel::PopulationPriorsSingle();
}

/**
 * Constructor to customize parameters.
 *
 * @param  population_mutation_rate Population mutation rate.
 * @param  germline_mutation_rate   Germline mutation rate.
 * @param  somatic_mutation_rate    Somatic mutation rate.
 * @param  sequencing_error_rate    Sequencing error rate.
 * @param  nucleotide_frequencies   Nucleotide frequencies changes the
 *                                  distribution of mutated nucleotides in the
 *                                  population priors.
 */
MultinomialTrioModel::MultinomialTrioModel(double population_mutation_rate,
                                           double germline_mutation_rate,
                                           double somatic_mutation_rate,
                                           double sequencing_error_rate,
                                           const RowVector4d &nucleotide_frequencies)
    : TrioModel(population_mutation_rate,
                germline_mutation_rate,
                somatic_mutation_rate,
                sequencing_error_rate,
                1.0,
                nucleotide_frequencies) {
  population_priors_ = TrioModel::PopulationPriors();
  population_priors_single_ = TrioModel::PopulationPriorsSingle();
}

/**
 * Returns 16 x 16 Eigen matrix. This is an order-relevant representation
 * of the possible events in the sample space that covers all possible parent
 * genotype combinations.
 *
 * Calls SpectrumProbability assuming infinite sites model using all enumerated
 * nucleotide counts at coverage 4x.
 *
 * @return  16 x 16 Eigen matrix in log e space where the (i, j) element is the
 *          probability that the mother has genotype i and the father has
 *          genotype j.
 */
Matrix16_16d MultinomialTrioModel::PopulationPriorsExpanded() {
  Matrix16_16d population_priors = Matrix16_16d::Zero();
  const Matrix16_16_4d kTwoParentCounts = TwoParentCounts();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      RowVector4d nucleotide_counts = kTwoParentCounts(i, j);
      double probability = SpectrumProbability(nucleotide_counts);
      population_priors(i, j) = probability;
    }
  }

  return population_priors;
}

/**
 * Returns the probability of the drawn alleles having a 4-0, 3-1, or 2-2
 * spectrum.
 *
 * @param  nucleotide_counts RowVector containing allele counts.
 * @return                   Probability of allele spectrum.
 */
double MultinomialTrioModel::SpectrumProbability(const RowVector4d &nucleotide_counts) {
  double p2_2 = population_mutation_rate_ / 2.0;
  double p3_1 = population_mutation_rate_ + population_mutation_rate_ / 3.0;
  double p4_0 = 1.0 - p3_1 - p2_2;

  if (IsInVector(nucleotide_counts, 4.0)) {
    return p4_0 * 0.25;  // p(4 allele) * p(position)
  } else if (IsInVector(nucleotide_counts, 3.0)) {
    return p3_1 * 0.25 / 3.0 * 0.25;  // p(3 allele) * p(1 allele) * p(position)
  } else if (IsInVector(nucleotide_counts, 2.0) &&
      !IsInVector(nucleotide_counts, 1.0)) {
    return p2_2 * 0.25 / 3.0 * 2.0 / 6.0;  // p(2 allele) * p(2 allele) * pair qualifier * p(position)
  } else {
    return 0.0;  // Counts do not match 4-0, 3-0, or 2-2 allele spectrum.
  }
}

/**
 * Calculates the probability of sequencing error for all read data. Assume
 * data contains 3 reads (child, mother, father). Assume the ReadDataVector is
 * already initialized in read_dependent_data_. Assume each chromosome is
 * equally likely to be sequenced.
 *
 * Adds the max element of all reads in ReadDataVector to
 * read_dependent_data_.max_elements before rescaling to normal space.
 */
void MultinomialTrioModel::SequencingProbabilityMat() {
  for (int read = 0; read < 3; ++read) {
    for (int genotype_idx = 0; genotype_idx < kGenotypeCount; ++genotype_idx) {
      auto alpha = alphas_.row(genotype_idx);
      // Converts alpha to double array.
      double p[kNucleotideCount] = {alpha(0), alpha(1), alpha(2), alpha(3)};
      // Converts read to unsigned int array.
      const ReadData &data = read_dependent_data_.read_data_vec[read];
      unsigned int n[kNucleotideCount] = {data.reads[0], data.reads[1],
                                          data.reads[2], data.reads[3]};
      double log_probability = gsl_ran_multinomial_lnpdf(kNucleotideCount, p, n);
      read_dependent_data_.sequencing_probability_mat(read, genotype_idx) = log_probability;
    }
  }
  
  // Rescales to normal space and records max element of all 3 reads together.
  double max_element = read_dependent_data_.sequencing_probability_mat.maxCoeff();
  read_dependent_data_.max_elements.push_back(max_element);

  // Calculates sequencing_probability_mat and splits into individual child
  // mother, and father vectors.
  read_dependent_data_.sequencing_probability_mat = exp(
    read_dependent_data_.sequencing_probability_mat.array() - max_element
  );
  read_dependent_data_.child_somatic_probability = read_dependent_data_.sequencing_probability_mat.row(0);
  read_dependent_data_.mother_somatic_probability = read_dependent_data_.sequencing_probability_mat.row(1);
  read_dependent_data_.father_somatic_probability = read_dependent_data_.sequencing_probability_mat.row(2);
}