/**
 * @file variant_visitor.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the VariantVisitor class.
 *
 * See top of variant_visitor.h for a complete description.
 */
#include "variant_visitor.h"


/**
 * Constructor calls base constructor PileupVisitor. base_cut cuts off
 * alignments for poor base quality. mapping_cut cuts off alignments for poor
 * mapping quality.
 *
 * @param  references Reference data from BAM file.
 * @param  header     SAM header data from BAM file.
 * @param  al         Alignment from BAM file.
 * @param  child_sm   Sample name of the child.
 * @param  mother_sm  Sample name of the mother.
 * @param  father_sm  Sample name of the father.
 */
VariantVisitor::VariantVisitor(const RefVector &references,
                               const string &child_sm,
                               const string &mother_sm,
                               const string &father_sm)
    : PileupVisitor(), references_{references}, child_sm_{child_sm},
      mother_sm_{mother_sm}, father_sm_{father_sm}, base_cut_{13},
      mapping_cut_{13} {
}

/**
 * Visits a BAM alignment at a certain position and counts the number of
 * nucleotide bases from read groups (RG) if it passes all cut values. The
 * sample name given through the command line is used here to identify which
 * sample it belongs to (child, mother, or father). Creates a ReadDataVector of
 * sequencing data and adds it to sites_.
 *
 * @param  pileupData  Pileup data.
 */
void VariantVisitor::Visit(const PileupPosition &pileupData) {
  ReadDataVector data_vec = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  string chr = references_[pileupData.RefId].RefName;

  for (auto it = begin(pileupData.PileupAlignments);
      it != end(pileupData.PileupAlignments); ++it) {
    if (it->Alignment.MapQuality >= mapping_cut_) {
      const int *pos = &it->PositionInAlignment;
      uint16_t bqual = static_cast<short>(it->Alignment.Qualities[*pos]) - 33;
      if (bqual >= base_cut_) {
        // Match individual identifier tag as child, mother, or father.
        string sm;
        it->Alignment.GetTag("RG", sm);
        int i = 0;

        if (sm.compare(child_sm_) == 0) {
          i = 0;  // Hard coded as implemented in TrioModel.
        } else if (sm.compare(mother_sm_) == 0) {
          i = 1;
        } else if (sm.compare(father_sm_) == 0) {
          i = 2;
        } else {
          cout << "ERROR: " << sm << " does not match." << endl;
          continue;
        }
        
        // Match nucleotide to individual at current position.
        char base = it->Alignment.QueryBases[*pos];
        uint16_t base_idx = ToNucleotideIndex(base);
        if (base_idx < 4) {
          data_vec[i].reads[base_idx]++;
        }
      }
    }
  }

  if (!HasZeroReadDataVector(data_vec)) {
    sites_.push_back(data_vec);
  }
}

TrioVector VariantVisitor::sites() const {
  return sites_;
}