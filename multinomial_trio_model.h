/**
 * @file multinomial_trio_model.h
 * @author Melissa Ip
 *
 * The MultinomialTrioModel class represents inherits TrioModel and uses
 * multinomial approximations and is based on the infinite sites model.
 */
#ifndef MULTINOMIALTRIOMODEL_H
#define MULTINOMIALTRIOMODEL_H

#include "trio_model.h"


class MultinomialTrioModel : public TrioModel {
 public:
  MultinomialTrioModel();  // Default constructor and constructor to customize parameters.
  MultinomialTrioModel(double population_mutation_rate,
                       double germline_mutation_rate,
                       double somatic_mutation_rate,
                       double sequencing_error_rate,
                       const RowVector4d &nucleotide_frequencies);

  double SpectrumProbability(const RowVector4d &nucleotide_counts);  // Calculates probability of allele spectrum given read counts.
  Matrix16_16d PopulationPriorsExpanded();
  void SequencingProbabilityMat();
};

#endif