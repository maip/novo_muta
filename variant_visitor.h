/**
 * @file variant_visitor.h
 * @author Melissa Ip
 *
 * The VariantVisitor class inherits the PileupVisitor class in BamTools.
 * It is used to parse BAM files and get sequencing reads as ReadDataVector
 * objects, which are accepted as input by the trio model.
 *
 * For further reference, see:
 * https://github.com/dwinter/accuMUlate
 */
#ifndef VARIANT_VISITOR_H
#define VARIANT_VISITOR_H

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "trio_model.h"

using namespace BamTools;


/**
 * VariantVisitor class header. See top of file for a complete description.
 */
class VariantVisitor : public PileupVisitor {
 public:
  VariantVisitor(const RefVector &references,
                 const string &child_sm,
                 const string &mother_sm,
                 const string &father_sm);
  ~VariantVisitor() {}
  void Visit(const PileupPosition &pileupData);
  TrioVector sites() const;

 private:
  RefVector references_;
  string child_sm_;
  string mother_sm_;
  string father_sm_;
  int base_cut_;
  int mapping_cut_;
  TrioVector sites_;
};

#endif
