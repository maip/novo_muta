/**
 * @file parameter_estimates.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the ParameterEstimates class.
 *
 * See top of parameter_estimates.h for a complete description.
 */
#include "parameter_estimates.h"


/**
 * Default constructor.
 *
 * @param  sites_count Number of sites.
 */
ParameterEstimates::ParameterEstimates(double sites_count)
    : theta_{0.0}, e_{0.0}, hom_{0.0}, het_{0.0}, som_{0.0}, germ_{0.0},
      n_s_{sites_count}, log_likelihood_{0.0},
      max_theta_{0.0}, max_e_{0.0}, count_{0} {
}

/**
 * Maximizes population mutation rate (theta). Calculated during the M-step of
 * expectation-maximization algorithm.
 */
double ParameterEstimates::MaxPopulationMutationRate() {
  return 6.0 / 11.0 * (1.0 - theta_/n_s_);
}

/**
 * Maximizes germline mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized germline mutation rate.
 */
double ParameterEstimates::MaxGermlineMutationRate() {
  double bracket_term = 1.0 - (4.0/3.0) * ((2*germ_) / n_s_);
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes somatic mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized somatic mutation rate.
 */
double ParameterEstimates::MaxSomaticMutationRate() {
  double bracket_term = 1.0 - (4.0/3.0) * ((6*som_) / n_s_);
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes sequencing error rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized sequencing error rate.
 */
double ParameterEstimates::MaxSequencingErrorRate() {
  double A = 2.0 * e_ + 2.0 * het_ + 2.0 * hom_;
  double B = -1.0 * (5.0 * e_ + 2.0 * het_ + 3.0 * hom_);
  double C = 3.0 * e_;
  double sqrt_term = sqrt(pow(B, 2) - 4.0*A*C);
  double first_term = -1.0 * B;
  double denominator = 2.0 * A;
  return (first_term - sqrt_term) / denominator;
}

/**
 * Calls the appropriate E-Step function for each summary statistic. Returns
 * true if successful, otherwise false.
 *
 * @param params  TrioModel object containing parameters.
 * @param sites   List of parsed trios.
 */
bool ParameterEstimates::Update(TrioModel &params, const TrioVector &sites) {
  Clear();  // Sets all statistics except number of sites to 0.
  for (const ReadDataVector data_vec : sites) {
    params.SetReadDependentData(data_vec);
    theta_ += SufficientStatistics::GetThetaStatistic(params);
    som_ += SufficientStatistics::GetSomaticStatistic(params);
    germ_ += SufficientStatistics::GetGermlineStatistic(params);
    e_ += SufficientStatistics::GetMismatchStatistic(params);
    hom_ += SufficientStatistics::GetHomozygousStatistic(params);
    het_ += SufficientStatistics::GetHeterozygousStatistic(params);

    log_likelihood_ += log(params.likelihood_read_dependent_data().denominator.sum);
    if (count_ == 0) {
      start_log_likelihood_ = log_likelihood_;
    }

    if (IsNan()) {
      cout << "ERROR: This site is nan." << endl;
      PrintReadDataVector(data_vec);
      Print(params.population_mutation_rate());
      // Print(params.sequencing_error_rate());
      Die("Cannot continue EM with a nan estimate.");
      return false;
    }
  }
  
  max_theta_ = MaxPopulationMutationRate();
  //max_e_ = MaxSequencingErrorRate();

  Print(params.population_mutation_rate());
  //Print(params.sequencing_error_rate());
  count_++;
  return true;
}

/**
 * Returns true if the sum of likelihood increases (and converges). To be
 * called at the end of the EM algorithm.
 *
 * @return  True if log likelihood increases for the data set, otherwise false.
 */
bool ParameterEstimates::IsLogLikelihoodIncreasing() {
  if (log_likelihood_ < start_log_likelihood_) {
    cout << "ERROR: Log likelihood is decreasing overall from "
         << start_log_likelihood_ << " to " << log_likelihood_ << endl;
    return false;
  }

  return true;
}

/**
 * Sets all statistics to 0 except n_s_, count_, and start_log_likelihood_.
 */
void ParameterEstimates::Clear() {
  theta_ = 0.0;
  e_ = 0.0;
  hom_ = 0.0;
  het_ = 0.0;
  som_ = 0.0;
  germ_ = 0.0;
  log_likelihood_ = 0.0;
  max_e_ = 0.0;
  max_theta_ = 0.0;
}

/**
 * Returns true if any of the statistics is NaN. Does not check for n_s_ or count_.
 */
bool ParameterEstimates::IsNan() {
  return (std::isnan(theta_) || std::isnan(e_) || std::isnan(hom_) ||
          std::isnan(het_) || std::isnan(som_) || std::isnan(germ_) ||
          std::isnan(log_likelihood_) || std::isnan(start_log_likelihood_)) ||
          std::isnan(max_e_) || std::isnan(max_theta_);
}

/**
 * Prints content.
 */
void ParameterEstimates::Print(double rate) {
  cout.precision(16);
  cout << "Iteration #" << count_ << " for " << n_s_ << " sites." << endl;
  cout << "^Theta:\t" << rate << endl
       // << "^E:\t"     << rate            << endl
       //<< "S_Som:\t"  << som_            << endl
       //<< "S_Germ:\t" << germ_           << endl
       // << "S_E:\t"    << e_              << endl
       // << "S_Hom:\t"  << hom_            << endl
       // << "S_Het:\t"  << het_            << endl
       << "Q_Log:\t"  << log_likelihood_ << endl
       << "~Theta:\t" << max_theta_ << endl << endl;
       // << "~E:\t"     << max_e_          << endl;
}

/**
 * Prints the maximum likelihood estimate of sequencing error rate after the
 * EM algorithm is done.
 */
void ParameterEstimates::PrintMaxPopulationMutationRateEstimate() {
  cout << "After " << count_ << " iterations of "
       << n_s_ << " sites:" << endl
       << "~Theta:\t" << max_theta_ << endl << endl;
}

/**
 * Prints the maximum likelihood estimate of sequencing error rate after the
 * EM algorithm is done.
 */
void ParameterEstimates::PrintMaxSequencingErrorRateEstimate() {
  cout << "After " << count_ << " iterations of "
       << n_s_ << " sites:" << endl
       << "~E:\t" << max_e_ << endl << endl;
}

double ParameterEstimates::theta() const {
  return theta_;
}

double ParameterEstimates::e() const {
  return e_;
}

double ParameterEstimates::hom() const {
  return hom_;
}

double ParameterEstimates::het() const {
  return het_;
}

double ParameterEstimates::som() const {
  return som_;
}

double ParameterEstimates::germ() const {
  return germ_;
}

double ParameterEstimates::n_s() const {
  return n_s_;
}

double ParameterEstimates::log_likelihood() const {
  return log_likelihood_;
}

double ParameterEstimates::max_theta() const {
  return max_theta_;
}

double ParameterEstimates::max_e() const {
  return max_e_;
}

int ParameterEstimates::count() const {
  return count_;
}

void ParameterEstimates::set_theta(double max) {
  theta_ = max;
}

void ParameterEstimates::set_e(double max) {
  e_ = max;
}

void ParameterEstimates::set_max_theta(double max) {
  max_theta_ = max;
}

void ParameterEstimates::set_max_e(double max) {
  max_e_ = max;
}

void ParameterEstimates::set_hom(double max) {
  hom_ = max;
}

void ParameterEstimates::set_het(double max) {
  het_ = max;
}

void ParameterEstimates::set_som(double max) {
  som_ = max;
}

void ParameterEstimates::set_germ(double max) {
  germ_ = max;
}

void ParameterEstimates::set_n_s(double num) {
  n_s_ = num;
}

void ParameterEstimates::set_log_likelihood(double num) {
  log_likelihood_ = num;
}
