/**
 * @file sufficient_statistics.h
 * @author Melissa Ip
 *
 * This file contains the implementation of the E-Step in the
 * expectation-maximization algorithm applied to the infinite sites model,
 * which will be updated as necessary in Fall 2014 and Spring 2015 to more
 * complex and biologically realistic models including the custom model using
 * Dirichlet-multinomial approximations instead of multinomial approximations.
 */
#ifndef SUFFICIENT_STATISTICS_H
#define SUFFICIENT_STATISTICS_H

#include "trio_model.h"


/**
 * SufficientStatistics namespace header. See top of file for a complete description.
 */
namespace SufficientStatistics {
  double GetPopulationMutationRateStatistic(const TrioModel &params);
  double GetSequencingErrorStatistic(const TrioModel &params,
                                     const Matrix3_16d &matches);
  double GetHomozygousStatistic(const TrioModel &params);
  double GetHeterozygousStatistic(const TrioModel &params);
  double GetMismatchStatistic(const TrioModel &params);
  RowVector16d GetHeterozygousMatches(const ReadData &data);
  Matrix3_16d GetHeterozygousMatches(const ReadDataVector &data_vec);
  RowVector16d GetHomozygousMatches(const ReadData &data);
  Matrix3_16d GetHomozygousMatches(const ReadDataVector &data_vec);
  RowVector16d GetMismatches(const ReadData &data);
  Matrix3_16d GetMismatches(const ReadDataVector &data_vec);
  double GetGermlineStatistic(const TrioModel &params);
  Matrix16_256d GermlineMutationCounts(const TrioModel &params);
  Matrix4_16d GermlineMutationCountsSingle(const TrioModel &params);
  double GetSomaticStatistic(const TrioModel &params);
  Matrix16_16d SomaticMutationCounts();
};

#endif
