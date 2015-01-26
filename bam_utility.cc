/**
 * @file bam_utility.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the BamUtility namespace.
 *
 * See top of bam_utility.h for a complete description.
 */
#include "bam_utility.h"


/**
 * Returns a list of trios parsed from input BAM file, matching reads given the
 * sample name of the child, mother, and father.
 *
 * @param  reader    BamReader instance.
 * @param  child_sm  Sample name of the child.
 * @param  mother_sm Sample name of the mother.
 * @param  father_sm Sample name of the father.
 * @return           List of trios parsed from input BAM file.
 */
TrioVector BamUtility::ParseSites(BamReader &reader,
                                  const string &child_sm,
                                  const string &mother_sm,
                                  const string &father_sm) {
  BamAlignment al;
  PileupEngine pileup;
  RefVector references = reader.GetReferenceData();
  VariantVisitor *v = new VariantVisitor(references, child_sm, mother_sm, father_sm);
  
  pileup.AddVisitor(v);
  while (reader.GetNextAlignment(al)) {
    pileup.AddAlignment(al);
  }

  pileup.Flush();
  return v->sites();
}

/**
 * Writes to text file the probability from each trio in sites.
 *
 * @param file_name Output file name.
 * @param sites     List of trios parsed from input BAM file.
 * @param params    TrioModel object containing parameters.
 */
void BamUtility::WriteProbability(const string &file_name,
                                  const TrioVector &sites,
                                  TrioModel &params) {
  vector<double> probabilities;
  for (auto site : sites) {
    double probability = params.MutationProbability(site);
    probabilities.push_back(probability);
  }

  ofstream fout(file_name);
  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();
}
