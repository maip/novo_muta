/**
 * @file pileup_driver.cc
 * @author Melissa Ip
 *
 * This file accepts 3 input pileup files. The pileup data is read and parsed
 * into sequencing read data that the TrioModel can process. Default parameter
 * values are used. Assume that all pileup files have the same number of sites
 * and thus can be aligned.
 * 
 * To compile on Herschel without using cmake and include GSL:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o pileup_driver pileup_utility.cc pileup_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./pileup_driver <output>.txt <child>.pileup <mother>.pileup <father>.pileup
 *
 * For example:
 * ./pileup_driver output.txt 878_21.pileup 892_21.pileup 891_21.pileup
 *
 * See top of pileup_utility.h for additional information.
 */
#include "em_algorithm.h" 
#include "pileup_utility.h"


int main(int argc, const char *argv[]) {
  if (argc < 5) {
    Die("USAGE: pileup_driver <output>.txt <child>.pileup <mother>.pileup "
        "<father>.pileup");
  }

  const string file_name = argv[1];
  const string child_pileup = argv[2];
  const string mother_pileup = argv[3];
  const string father_pileup = argv[4];

  ifstream child(child_pileup);
  ifstream mother(mother_pileup);
  ifstream father(father_pileup);
  if (!child.is_open() || 0 != child.fail() || !mother.is_open() ||
      0 != mother.fail() || !father.is_open() || 0 != father.fail()) {
    Die("Input file cannot be read.");
  }  

  TrioVector sites = PileupUtility::ParseSites(child, mother, father);
  
  // Feature: Write probability to output.
  // vector<double> probabilities = PileupUtility::GetProbability(sites);
  // PileupUtility::WriteProbability(file_name, probabilities);

  // Feature: EM algorithm.
  cout.precision(16);
  TrioModel params;
  ParameterEstimates *stats = EstimateParameters(params, sites);
  stats->IsLogLikelihoodIncreasing();
  stats->PrintMaxSequencingErrorRateEstimate();

  child.close();
  mother.close();
  father.close();

  return 0;
}
