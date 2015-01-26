/**
 * @file bam_utility.h
 * @author Melissa Ip
 *
 * This file contains methods to parse sites from BAM files into a TrioVector
 * object and to write the probabilities of the trios to a text file.
 */
#ifndef BAM_UTILITY_H
#define BAM_UTILITY_H

#include <fstream>

#include "variant_visitor.h"


/**
 * BamUtility namespace header. See top of file for a complete description.
 */
namespace BamUtility {
  TrioVector ParseSites(BamReader &reader, const string &child_sm,
                        const string &mother_sm, const string &father_sm);
  void WriteProbability(const string &file_name, const TrioVector &sites,
                        TrioModel &params);
};

#endif
