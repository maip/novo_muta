/**
 * @file trio_model.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the TrioModel class.
 *
 * See top of trio_model.h for a complete description.
 */
#include "trio_model.h"


/**
 * Default constructor.
 * 
 * sequencing_probability_mat is created or updated if sequencing_error_rate_
 * or dirichlet_dispersion_ is changed when MutationProbability() or
 * SetReadDependentData() is called. dirichlet_dispersion_ is not used in the
 * infinite sites model version.
 */
TrioModel::TrioModel()
    : population_mutation_rate_{0.001},
      germline_mutation_rate_{2e-8},
      somatic_mutation_rate_{2e-8},
      sequencing_error_rate_{0.005},
      dirichlet_dispersion_{1000.0},
      nucleotide_frequencies_{0.25, 0.25, 0.25, 0.25} {
  population_priors_ = PopulationPriors();
  population_priors_single_ = PopulationPriorsSingle();
  SetGermlineMutationProbabilities();
  germline_probability_mat_single_ = GermlineProbabilityMatSingle();
  germline_probability_mat_ = GermlineProbabilityMat();
  germline_probability_mat_num_ = GermlineProbabilityMat(true);
  somatic_probability_mat_ = SomaticProbabilityMat();
  somatic_probability_mat_diag_ = SomaticProbabilityMatDiag();
  alphas_ = GetAlphas();
}

/**
 * Constructor to customize parameters.
 *
 * @param  population_mutation_rate Population mutation rate.
 * @param  germline_mutation_rate   Germline mutation rate.
 * @param  somatic_mutation_rate    Somatic mutation rate.
 * @param  sequencing_error_rate    Sequencing error rate.
 * @param  dirichlet_dispersion     Dirichlet dispersion changes the magnitude
 *                                  of the probability of mutation.
 * @param  nucleotide_frequencies   Nucleotide frequencies changes the
 *                                  distribution of mutated nucleotides in the
 *                                  population priors.
 */
TrioModel::TrioModel(double population_mutation_rate,
                     double germline_mutation_rate,
                     double somatic_mutation_rate,
                     double sequencing_error_rate,
                     double dirichlet_dispersion,
                     const RowVector4d &nucleotide_frequencies)
    : population_mutation_rate_{population_mutation_rate},
      germline_mutation_rate_{germline_mutation_rate},
      somatic_mutation_rate_{somatic_mutation_rate},
      sequencing_error_rate_{sequencing_error_rate},
      dirichlet_dispersion_{dirichlet_dispersion},
      nucleotide_frequencies_{nucleotide_frequencies} {
  population_priors_ = PopulationPriors();
  population_priors_single_ = PopulationPriorsSingle();
  SetGermlineMutationProbabilities();
  germline_probability_mat_single_ = GermlineProbabilityMatSingle();
  germline_probability_mat_ = GermlineProbabilityMat();
  germline_probability_mat_num_ = GermlineProbabilityMat(true);
  somatic_probability_mat_ = SomaticProbabilityMat();
  somatic_probability_mat_diag_ = SomaticProbabilityMatDiag();
  alphas_ = GetAlphas();
}

/**
 * Implements the trio model for a single site by calling the functions on the
 * left in the following diagram. The function names label the arrow-denoted
 * processes in the population model.
 *
 *                             Population          Population
 * PopulationPriors()              |                   |
 *                                 v                   v
 *                               Mother              Father
 *                               Zygotic             Zygotic
 *                               Diploid             Diploid
 *                               Genotype            Genotype
 * GermlineMutation()              |   \             / |  
 * GermlineProbabilityMat()        |    v           v  |
 *                                 |      Daughter     |
 *                                 |      Zygotic      |
 *                                 |      Diploid      |
 *                                 |      Genotype     |
 * SomaticMutation()               |         |         |
 * SomaticProbabilityMat()         v         v         v
 *                               Mother   Daughter  Father
 *                               Somatic  Somatic   Somatic
 *                               Diploid  Diploid   Diploid
 *                               Genotype Genotype  Genotype
 * SequencingProbabilityMat()      |         |         |
 *                                 v         v         v
 *                               Mother   Daughter  Father
 *                               Genotype Genotype  Genotype
 *                               Reads    Reads     Reads
 *
 * ReadDependentData is initialized in SetReadDependentData() call.
 *
 * @param   data_vec Read counts in order of child, mother and father.
 * @return           Probability of mutation given read data and parameters.
 */
double TrioModel::MutationProbability(const ReadDataVector &data_vec) {
  SetReadDependentData(data_vec);

  return 1 - (read_dependent_data_.numerator.sum /
              read_dependent_data_.denominator.sum);
}

/**
 * Initializes and updates read_dependent_data_.sequencing_probability_mat and
 * individual somatic probabilities using sequencing_error_rate_ as well as the
 * other probabilities along the branch.
 *
 * Follows the model diagram in MutationProbability.
 *
 * @param   data_vec Read counts in order of child, mother and father.
 */
void TrioModel::SetReadDependentData(const ReadDataVector &data_vec) {
  read_dependent_data_ = ReadDependentData(data_vec);  // First intialized.

  SequencingProbabilityMat();
  SomaticTransition();
  GermlineTransition();
  SomaticTransition(true);
  GermlineTransition(true);
}

/**
 * Returns 1 x 256 Eigen probability RowVector. This is an order-relevant
 * representation of the possible events in the sample space that covers all
 * possible parent genotype combinations. For example:
 *
 * [P(AAAA), P(AAAC), P(AAAG), P(AAAT), P(AACA), P(AACC), P(AACG)...]
 *
 * Resizes the original 16 x 16 matrix to 1 x 256.
 *
 * @return  1 x 256 Eigen probability RowVector in log e space where the i
 *          element is a unique parent pair genotype.
 */
RowVector256d TrioModel::PopulationPriors() {
  RowVector256d population_priors_flattened;
  Matrix16_16d population_priors_expanded = PopulationPriorsExpanded();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      int idx = i * kGenotypeCount + j;
      population_priors_flattened(idx) = population_priors_expanded(i, j);
    }
  }

  return population_priors_flattened;
}

/**
 * Returns 16 x 16 Eigen matrix. This is an order-relevant representation
 * of the possible events in the sample space that covers all possible parent
 * genotype combinations.
 *
 * Creates nucleotide mutation frequencies {alpha_A, alpha_C, alpha_G, alpha_T}
 * based on the nucleotide frequencies and population mutation rate (theta).
 * These frequencies and nucleotide counts {n_A, n_C, n_G, n_T} are used in the
 * Dirichlet multinomial.
 *
 * For example, both parents have genotype AA, resulting in N = 4:
 *
 * population_mutation_rate_ = 0.00025;
 * nucleotide_frequencies_ << 0.25, 0.25, 0.25, 0.25;
 * nucleotide_counts = {4, 0, 0, 0};
 *
 * @return  16 x 16 Eigen matrix in log e space where the (i, j) element is the
 *          probability that the mother has genotype i and the father has
 *          genotype j.
 */
Matrix16_16d TrioModel::PopulationPriorsExpanded() {
  // Calculates nucleotide mutation frequencies using given mutation rate.
  RowVector4d nucleotide_mutation_frequencies = (nucleotide_frequencies_ *
    population_mutation_rate_);
  Matrix16_16d population_priors = Matrix16_16d::Zero();
  const Matrix16_16_4d kTwoParentCounts = TwoParentCounts();

  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      // Convert nucleotide_counts to ReadData for DirichletMultinomialLog().
      RowVector4d nucleotide_counts = kTwoParentCounts(i, j);
      ReadData nucleotide_read;
      for (int k = 0; k < kNucleotideCount; ++k) {
        nucleotide_read.reads[k] = nucleotide_counts(k);
      }
      // Calculates probability using the Dirichlet multinomial in normal space.
      double log_probability = DirichletMultinomialLog(
        nucleotide_mutation_frequencies,
        nucleotide_read
      );
      population_priors(i, j) = exp(log_probability);
    }
  }

  return population_priors;
}

/**
 * Returns 1 x 16 Eigen RowVector population priors for a single parent.
 */
RowVector16d TrioModel::PopulationPriorsSingle() {
  return PopulationPriorsExpanded().rowwise().sum();
}

/**
 * Calculates set of possible germline mutation probabilities given germline
 * mutation rate. Weighted based on if parent genotype is homozygous or
 * heterozygous.
 */
void TrioModel::SetGermlineMutationProbabilities() {
  double exp_term = exp(-4.0/3.0 * germline_mutation_rate_);
  homozygous_match_ = 0.25 + 0.75 * exp_term;
  heterozygous_match_ = 0.25 + 0.25 * exp_term;
  mismatch_ = 0.25 - 0.25 * exp_term;
}

/**
 * Calculates the probability of germline mutation and parent chromosome
 * donation by default. Assume the first chromosome is associated with
 * the mother and the second chromosome is associated with the father.
 *
 * @param  child_nucleotide_idx Index of child allele.
 * @param  parent_genotype_idx  Index of parent genotype.
 * @param  no_mutation_flag     False by default. Set to true to calculate
 *                              probability of no mutation.
 * @return                      Probability of germline mutation.
 */
double TrioModel::GermlineMutation(int child_nucleotide_idx,
                                   int parent_genotype_idx,
                                   bool no_mutation_flag) {
  // Determines if the comparison is homozygous, heterozygous or no match.
  if (IsAlleleInParentGenotype(child_nucleotide_idx, parent_genotype_idx)) {
    if (parent_genotype_idx % 5 == 0) {  // Homozygous genotypes are divisible by 5.
      return homozygous_match_;
    } else {
      if (no_mutation_flag) {
        return homozygous_match_ / 2;
      } else {
        return heterozygous_match_;
      }
    }
  } else {
    if (no_mutation_flag) {
      return 0.0;
    } else {
      return mismatch_;
    }
  }
}

/**
 * Calculates the germline probability matrix for the offspring using the given
 * mutation rate, derived from the Kronecker product of a 4 x 16 Eigen 
 * germline probability matrix for a single parent with itself.
 *
 * This is a transition matrix used to mutate the child germline genotype.
 *
 * @param  no_mutation_flag False by default. Set to true to calculate
 *                          probability of no mutation.
 * @return                  16 x 256 Eigen probability matrix.
 */
Matrix16_256d TrioModel::GermlineProbabilityMat(bool no_mutation_flag) {
  return KroneckerProduct(GermlineProbabilityMatSingle(no_mutation_flag));
}

/**
 * Calculates the germline probability matrix for the offspring using the given
 * mutation rate.
 *
 * This is a transition matrix used to mutate the child germline genotype, that
 * is the probability of one random allele from one parent being mutated in the
 * germline. The current trio model will use the Kronecker product of this
 * matrix to produce the transition matrix of germline probabilities accounting
 * for both parents.
 *
 * @param  no_mutation_flag False by default. Set to true to calculate
 *                          probability of no mutation.
 * @return                  4 x 16 Eigen probability matrix.
 */
Matrix4_16d TrioModel::GermlineProbabilityMatSingle(bool no_mutation_flag) {
  Matrix4_16d germline_probability_mat = Matrix4_16d::Zero();
  for (int i = 0; i < kNucleotideCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      double probability = GermlineMutation(i, j, no_mutation_flag);
      germline_probability_mat(i, j) = probability;
    }
  }

  return germline_probability_mat;
}

/**
 * Calculates the probability of somatic mutation.
 *
 * @param  nucleotide_idx       Index of nucleotide.
 * @param  other_nucleotide_idx Index of another nucleotide to be compared.
 * @return                      Probability of somatic mutation.
 */
double TrioModel::SomaticMutation(int nucleotide_idx, int other_nucleotide_idx) {
  double exp_term = exp(-4.0/3.0 * somatic_mutation_rate_);
  double term = 0.25 * (1 - exp_term);
  if (nucleotide_idx == other_nucleotide_idx) {  // Indicator function.
    return term + exp_term;
  } else {
    return term;
  }
}

/**
 * Computes event space for somatic nucleotide given a genotype nucleotide for
 * a single chromosome. Combines event spaces for two chromosomes independent
 * of each other and calculates somatic mutation probability matrix for a 
 * single parent.
 *
 * This is a transition matrix used to mutate the somatic genotypes.
 *
 * @return  16 x 16 Eigen probability matrix where the first dimension is the
 *          original somatic genotypes and the second dimension is the mutated 
 *          genotype.
 */
Matrix16_16d TrioModel::SomaticProbabilityMat() {
  Matrix4d somatic_probability_mat = Matrix4d::Zero();
  for (int i = 0; i < kNucleotideCount; ++i) {
    for (int j = 0; j < kNucleotideCount; ++j) {
      double probability = SomaticMutation(i, j);
      somatic_probability_mat(i, j) = probability;
    }
  }
  return KroneckerProduct(somatic_probability_mat);
}

/**
 * Returns diagonal of somatic_probability_mat_. Assumes SomaticProbabilityMat
 * has been called.
 *
 * @return  16 x 16 Eigen matrix diagonal of somatic_probability_mat_.
 */
Matrix16_16d TrioModel::SomaticProbabilityMatDiag() {
  return somatic_probability_mat_.diagonal().asDiagonal();
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
void TrioModel::SequencingProbabilityMat() {
  for (int read = 0; read < 3; ++read) {
    for (int genotype_idx = 0; genotype_idx < kGenotypeCount; ++genotype_idx) {
      auto alpha = alphas_.row(genotype_idx);
      const ReadData &data = read_dependent_data_.read_data_vec[read];
      double log_probability = DirichletMultinomialLog(alpha, data);
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

/**
 * Multiplies sequencing probability vectors by somatic transition matrix.
 *
 * @param  is_numerator True if calculating probability of numerator.
 */
void TrioModel::SomaticTransition(bool is_numerator) {
  if (!is_numerator) {
    read_dependent_data_.denominator.child_zygotic_probability = (
      read_dependent_data_.child_somatic_probability * somatic_probability_mat_
    );
    read_dependent_data_.denominator.mother_zygotic_probability = (
      read_dependent_data_.mother_somatic_probability * somatic_probability_mat_
    );
    read_dependent_data_.denominator.father_zygotic_probability = (
      read_dependent_data_.father_somatic_probability * somatic_probability_mat_
    );
  } else {
    read_dependent_data_.numerator.child_zygotic_probability = (
      read_dependent_data_.child_somatic_probability * somatic_probability_mat_diag_
    );
    read_dependent_data_.numerator.mother_zygotic_probability = (
      read_dependent_data_.mother_somatic_probability * somatic_probability_mat_diag_
    );
    read_dependent_data_.numerator.father_zygotic_probability = (
      read_dependent_data_.father_somatic_probability * somatic_probability_mat_diag_
    );
  }
}

/**
 * Calculates denominator, probability of the observed data or numerator,
 * probability of no mutation.
 *
 * @param  is_numerator True if calculating probability of numerator.
 */
void TrioModel::GermlineTransition(bool is_numerator) {
  if (!is_numerator) {
    read_dependent_data_.denominator.child_germline_probability = (
      read_dependent_data_.denominator.child_zygotic_probability *
      germline_probability_mat_
    );
    read_dependent_data_.denominator.parent_probability = KroneckerProduct(
      read_dependent_data_.denominator.mother_zygotic_probability,
      read_dependent_data_.denominator.father_zygotic_probability
    );
    read_dependent_data_.denominator.root_mat = GetRootMat(
      read_dependent_data_.denominator.child_germline_probability,
      read_dependent_data_.denominator.parent_probability
    );
    read_dependent_data_.denominator.sum = read_dependent_data_.denominator.root_mat.sum();
  } else {
    read_dependent_data_.numerator.child_germline_probability = (
      read_dependent_data_.numerator.child_zygotic_probability *
      germline_probability_mat_num_
    );
    read_dependent_data_.numerator.parent_probability = KroneckerProduct(
      read_dependent_data_.numerator.mother_zygotic_probability,
      read_dependent_data_.numerator.father_zygotic_probability
    );
    read_dependent_data_.numerator.root_mat = GetRootMat(
      read_dependent_data_.numerator.child_germline_probability,
      read_dependent_data_.numerator.parent_probability
    );
    read_dependent_data_.numerator.sum = read_dependent_data_.numerator.root_mat.sum();
  }
}

/**
 * Returns root matrix for use in MutationProbability.
 *
 * @param  child_germline_probability Matrix after child germline transition.
 * @param  parent_probability         Kronecker product of both parent vectors.
 * @return                            1 x 256 final matrix at the root of tree.
 */
RowVector256d TrioModel::GetRootMat(const RowVector256d &child_germline_probability,
                                    const RowVector256d &parent_probability) {
  return child_germline_probability.cwiseProduct(
    parent_probability).cwiseProduct(population_priors_);
}

/**
 * Generates a 16 x 4 alpha frequencies matrix given the sequencing error rate
 * and dirichlet dispersion. The order of the alpha frequencies correspond to
 * the genotypes. Each alpha should sum to 1.
 * 
 * Current values are placeholders until they are estimated in Spring 2014.
 *
 * @return  16 x 4 Eigen matrix of Dirichlet multinomial alpha parameters
 *          alpha = (alpha_1, ..., alpha_K) for a K-category Dirichlet
 *          distribution (where K = 4 = kNucleotideCount) that vary with each
 *          combination of parental genotype and reference nucleotide.
 */
Matrix16_4d TrioModel::GetAlphas() {
  Matrix16_4d alphas;
  double homozygous = 1.0 - sequencing_error_rate_;
  double mismatch = sequencing_error_rate_ / 3.0;
  double heterozygous = 0.5 - mismatch;

  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kNucleotideCount; ++j) {
      if (IsAlleleInParentGenotype(j, i)) {
        if (i % 5 == 0) {  // Homozygous genotypes are divisible by 5.
          alphas(i, j) = homozygous;
        } else {
          alphas(i, j) = heterozygous;
        }
      } else {
        alphas(i, j) = mismatch;
      }
    }
  }

  return alphas * dirichlet_dispersion_;
}

/**
 * Returns true if the two TrioModel objects are equal to each other within
 * epsilon precision.
 *
 * If the underflow/overflow bug is fixed, it would be more appropriate to
 * remove the approximation comparison and check for direct equality.
 * 
 * @param  other TrioModel object to be compared.
 * @return       True if the two TrioModel objects are equal to each other.
 */
bool TrioModel::Equals(const TrioModel &other) {
  bool attr_table[12] = {
    Equal(population_mutation_rate_, other.population_mutation_rate_),
    Equal(germline_mutation_rate_, other.germline_mutation_rate_),
    Equal(somatic_mutation_rate_, other.germline_mutation_rate_),
    Equal(sequencing_error_rate_, other.sequencing_error_rate_),
    Equal(dirichlet_dispersion_, other.dirichlet_dispersion_),
    nucleotide_frequencies_.isApprox(other.nucleotide_frequencies_, kEpsilon),
    population_priors_single_.isApprox(other.population_priors_single_, kEpsilon),
    population_priors_.isApprox(other.population_priors_, kEpsilon),
    germline_probability_mat_single_.isApprox(
        other.germline_probability_mat_single_,
        kEpsilon),
    germline_probability_mat_.isApprox(
        other.germline_probability_mat_,
        kEpsilon),
    somatic_probability_mat_.isApprox(other.somatic_probability_mat_, kEpsilon),
    read_dependent_data_.sequencing_probability_mat.isApprox(
        other.read_dependent_data_.sequencing_probability_mat,
        kEpsilon)
  };

  if (all_of(begin(attr_table), end(attr_table), [](bool i) { return i; })) {
    return true;
  } else {
    return false;
  }
}

double TrioModel::population_mutation_rate() const {
  return population_mutation_rate_;
}

/**
 * Sets population_mutation_rate_, population_priors_single_ and population_priors_.
 */
void TrioModel::set_population_mutation_rate(double rate) {
  population_mutation_rate_ = rate;
  population_priors_ = PopulationPriors();
  population_priors_single_ = PopulationPriorsSingle();
}

double TrioModel::germline_mutation_rate() const {
  return germline_mutation_rate_;
}

/**
 * Sets germline_mutation_rate_, germline_probability_mat_single,
 * germline_probability_mat_ and germline_probability_mat_num_.
 */
void TrioModel::set_germline_mutation_rate(double rate) {
  germline_mutation_rate_ = rate;
  SetGermlineMutationProbabilities();
  germline_probability_mat_single_ = GermlineProbabilityMatSingle();
  germline_probability_mat_ = GermlineProbabilityMat();
  germline_probability_mat_num_ = GermlineProbabilityMat(true);
}

double TrioModel::homozygous_match() const {
  return homozygous_match_;
}

double TrioModel::heterozygous_match() const {
  return heterozygous_match_;
}

double TrioModel::mismatch() const {
  return mismatch_;
}

double TrioModel::somatic_mutation_rate() const {
  return somatic_mutation_rate_;
}

/**
 * Sets somatic_mutation_rate_, somatic_probability_mat_ and
 * somatic_probability_mat_diag_.
 */
void TrioModel::set_somatic_mutation_rate(double rate) {
  somatic_mutation_rate_ = rate;
  somatic_probability_mat_ = SomaticProbabilityMat();
  somatic_probability_mat_diag_ = SomaticProbabilityMatDiag();
}

double TrioModel::sequencing_error_rate() const {
  return sequencing_error_rate_;
}

/**
 * Sets sequencing_error_rate_ and alphas_.
 */
void TrioModel::set_sequencing_error_rate(double rate) {
  sequencing_error_rate_ = rate;
  alphas_ = GetAlphas();
}

double TrioModel::dirichlet_dispersion() const {
  return dirichlet_dispersion_;
}

/**
 * Sets dirichlet_dispersion_ and alphas_.
 */
void TrioModel::set_dirichlet_dispersion(double dispersion) {
  dirichlet_dispersion_ = dispersion;
  alphas_ = GetAlphas();
}

RowVector4d TrioModel::nucleotide_frequencies() const {
  return nucleotide_frequencies_;
}

/**
 * Sets nucleotide_frequencies_, population_priors_ and population_priors_single_.
 */
void TrioModel::set_nucleotide_frequencies(const RowVector4d &frequencies) {
  nucleotide_frequencies_ = frequencies;
  population_priors_ = PopulationPriors();
  population_priors_single_ = PopulationPriorsSingle();
}

RowVector16d TrioModel::population_priors_single() const {
  return population_priors_single_;
}

RowVector256d TrioModel::population_priors() const {
  return population_priors_;
}

Matrix4_16d TrioModel::germline_probability_mat_single() const {
  return germline_probability_mat_single_;
}

Matrix16_256d TrioModel::germline_probability_mat() const {
  return germline_probability_mat_;
}

Matrix16_256d TrioModel::germline_probability_mat_num() const {
  return germline_probability_mat_num_;
}

Matrix16_16d TrioModel::somatic_probability_mat() const {
  return somatic_probability_mat_;
}

Matrix16_16d TrioModel::somatic_probability_mat_diag() const {
  return somatic_probability_mat_diag_;
}

Matrix3_16d TrioModel::sequencing_probability_mat() const {
  return read_dependent_data_.sequencing_probability_mat;
}

Matrix16_4d TrioModel::alphas() const {
  return alphas_;
}

ReadDependentData TrioModel::read_dependent_data() const {
  return read_dependent_data_;
}
