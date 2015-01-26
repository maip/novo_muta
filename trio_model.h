/**
 * @file trio_model.h
 * @author Melissa Ip
 *
 * The TrioModel class is used to create an object that contains sequencing
 * read data, estimation parameters, and useful matrices in calculating the
 * probability of de novo mutation.
 *
 * The probability of mutation is calculated using a modified trio model,
 * similar but different to the model described in the following paper:
 * 
 * Cartwright et al.: Family-Based Method for Capturing De Novo Mutations
 * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728889/
 *
 * This TrioModel class uses Dirichlet-multinomial approximations while the
 * derived class MultinomialTrioModel uses simpler multinomial approximations
 * and is based on the infinite sites model.
 *
 * Example usage:
 *
 *   TrioModel params;  // Uses default parameters with Dirichlet-multinomial.
 *   ReadDataVector data = {  // Sequencing data in order: child, mother, father.
 *     {30, 0, 0, 0},
 *     {30, 0, 0, 0},
 *     {30, 0, 0, 0}
 *   };
 *
 *   double probability = params.MutationProbability(data);
 *   data.set_germline_mutation_rate(0.000001);
 *   double new_probability = params.MutationProbability(data);
 */
#ifndef TRIO_MODEL_H
#define TRIO_MODEL_H

#include "read_dependent_data.h"


/**
 * TrioModel class header. See top of file for a complete description.
 */
class TrioModel {
 public:
  TrioModel();  // Default constructor and constructor to customize parameters.
  TrioModel(double population_mutation_rate,
            double germline_mutation_rate,
            double somatic_mutation_rate,
            double sequencing_error_rate,
            double dirichlet_dispersion,
            const RowVector4d &nucleotide_frequencies);
  double MutationProbability(const ReadDataVector &data_vec);  // Calculates probability of mutation given input read data.
  void SetReadDependentData(const ReadDataVector &data_vec);
  bool Equals(const TrioModel &other);  // True if the two TrioModel objects are equal to each other.
  double population_mutation_rate() const;  // Get and set functions.
  void set_population_mutation_rate(double rate);
  double germline_mutation_rate() const;
  void set_germline_mutation_rate(double rate);
  double homozygous_match() const;
  double heterozygous_match() const;
  double mismatch() const;
  double somatic_mutation_rate() const;
  void set_somatic_mutation_rate(double rate);
  double sequencing_error_rate() const;
  void set_sequencing_error_rate(double rate);
  double dirichlet_dispersion() const;
  void set_dirichlet_dispersion(double dispersion);
  RowVector4d nucleotide_frequencies() const;
  void set_nucleotide_frequencies(const RowVector4d &frequencies);
  RowVector16d population_priors_single() const;
  RowVector256d population_priors() const;
  Matrix4_16d germline_probability_mat_single() const;
  Matrix16_256d germline_probability_mat() const;
  Matrix16_256d germline_probability_mat_num() const;
  Matrix16_16d somatic_probability_mat() const;
  Matrix16_16d somatic_probability_mat_diag() const;
  Matrix3_16d sequencing_probability_mat() const;
  Matrix16_4d alphas() const;
  ReadDependentData read_dependent_data() const;

 protected:
  void GermlineTransition(bool is_numerator=false);  // Helper functions for MutationProbability.
  void SomaticTransition(bool is_numerator=false);
  RowVector256d GetRootMat(const RowVector256d &child_germline_probability,
                           const RowVector256d &parent_probability);
  RowVector256d PopulationPriors();  // Functions for setting up the model and relevant arrays.
  virtual Matrix16_16d PopulationPriorsExpanded();
  RowVector16d PopulationPriorsSingle();
  void SetGermlineMutationProbabilities();
  double GermlineMutation(int child_nucleotide_idx, int parent_genotype_idx,
                          bool no_mutation_flag);
  Matrix4_16d GermlineProbabilityMatSingle(bool no_mutation_flag=false);
  Matrix16_256d GermlineProbabilityMat(bool no_mutation_flag=false);
  double SomaticMutation(int nucleotide_idx, int other_nucleotide_idx);
  Matrix16_16d SomaticProbabilityMat();
  Matrix16_16d SomaticProbabilityMatDiag();
  virtual void SequencingProbabilityMat();
  Matrix16_4d GetAlphas();

  // Instance member variables.
  double population_mutation_rate_;
  double homozygous_match_;
  double heterozygous_match_;
  double mismatch_;
  double germline_mutation_rate_;
  double somatic_mutation_rate_;
  double sequencing_error_rate_;
  double dirichlet_dispersion_;  // Set to 1.0 for MultinomialTrioModel.
  RowVector4d nucleotide_frequencies_;
  Matrix16_4d alphas_;
  RowVector16d population_priors_single_;  // Unused.
  RowVector256d population_priors_;
  Matrix4_16d germline_probability_mat_single_;
  Matrix16_256d germline_probability_mat_;
  Matrix16_256d germline_probability_mat_num_;
  Matrix16_16d somatic_probability_mat_;
  Matrix16_16d somatic_probability_mat_diag_;
  ReadDependentData read_dependent_data_;  // Contains TreePeel class.
};

#endif
