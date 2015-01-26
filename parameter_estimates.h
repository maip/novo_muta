/**
 * @file parameter_estimates.h
 * @author Melissa Ip
 *
 * The ParameterEstimates class is used to hold ~S_T of type T sufficient
 * statistics produced in the E-Step of the expectation-maximization algorithm.
 * It also contains the sum of each sufficient statistic after each
 * site/iteration. This calculates the maximized-likelihood estimates for the
 * M-Step.
 */
#ifndef PARAMETER_ESTIMATES_H
#define PARAMETER_ESTIMATES_H

#include "sufficient_statistics.h"


/**
 * ParameterEstimates class header. See top of file for a complete description.
 */
class ParameterEstimates {
 public:
  ParameterEstimates(double sites_count);  // Default constructor.
  ~ParameterEstimates() {}
  bool Update(TrioModel &params, const TrioVector &sites);
  void Clear();
  void Print();

  double MaxGermlineMutationRate();  // M-step functions.
  double MaxSomaticMutationRate();
  double MaxSequencingErrorRate();
  bool IsLogLikelihoodIncreasing();
  void PrintMaxSequencingErrorRateEstimate();
  bool IsNan();

  double e() const;  // Get and set functions.
  double hom() const;
  double het() const;
  double som() const;
  double germ() const;
  double n_s() const;
  double log_likelihood() const;
  double start_log_likelihood() const;
  double max_e() const;
  int count() const;
  void set_e(double max);
  void set_hom(double max);
  void set_het(double max);
  void set_som(double max);
  void set_germ(double max);
  void set_n_s(double num);
  void set_log_likelihood(double num);

 private:
  // Instance member variables.
  double e_;                     // Expected statistic for sequencing error rate.
  double hom_;                   // Expected statistic for homozygous matches.
  double het_;                   // Expected statistic for heterozygous matches.
  double som_;                   // Expected statistic for somatic mutation.
  double germ_;                  // S_M + S_F. Expected statistic for germline mutation.
  double n_s_;                   // Number of sites.
  double log_likelihood_;        // Log likelihood of P(R,H)/denominator sum.
  double start_log_likelihood_;  // First calculated log likelihood.
  double max_e_;                 // Maximum likehood estimate for sequencing error rate.
  int count_;                    // Number of iterations in EM algorithm performed.
};

#endif
