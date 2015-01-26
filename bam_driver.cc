/**
 * @file bam_driver.cc
 * @author Melissa Ip
 *
 * This is the driver for parsing a BAM file. The input file must contain
 * all reads for the trio (child, mother, father) and have the appropriate tags.
 *
 * SM refers to the name of the sample that identifies it as belonging to the
 * child, mother, or father. This will vary depending on the data and must be
 * known in order to parse out the read data.
 * 
 * Useful commands to view overview, header, and a certain region:
 * samtools idxstats <name>.bam
 * samtools view -H <name>.bam <chr>:<pos1>-<pos2>
 * samtools view <name>.bam <chr>:<pos1>-<pos2>
 * samtools tview -p <chr>:<pos1>-<pos2> <name>.bam
 *
 * To slice a region out and merge it with other slices:
 * samtools view -b <name>.bam <chr>:<pos1>-<pos2> > <output>.bam
 *
 * To output all header data that belong to a certain region:
 * samtools view -o /home/mip/<header>.txt -H <name>.bam <chr>:<pos1>-<pos2>
 *
 * Make read group file for bam1 to bamn, tab separated. <header>.txt can be
 * parsed by filtering out lines that begin with @RG with parse_bam_header.cc.
 *
 * Skip this step if there is only one bam file that contains the bam
 * alignments for child, mother, and father.
 * To merge separate bam files for each family member:
 * samtools merge -rh <rg>.txt <output>.bam <bam1>.bam <bamn>.bam
 *
 * To prepare the output for the index:
 * samtools sort <output>.bam <output_sorted>.bam
 *
 * To create the index:
 * samtools index <output_sorted>.bam <output>.index
 *
 * <output_sorted>.bam is accepted as input for a specified region.
 * 
 * To compile on Herschel without using cmake and include GSL and BamTools:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/mip/novo_muta_infinite_sites_model/bamtools/lib -I/home/mip/novo_muta_infinite_sites_model/bamtools/include -lbamtools -I/home/mip/novo_muta_infinite_sites_model/bamtools/src -lbamtools-utils -I/usr/local/include -o bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc sufficient_statistics.cc parameter_estimates.cc em_algorithm.cc bam_utilities.cc bam_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./bam_driver <input>.bam <child SM> <mother SM> <father SM>
 *
 * For example:
 * ./bam_driver output_sorted.bam NA12828 NA12892 NA12891
 */
#include "bam_utility.h"
#include "em_algorithm.cc"


int main(int argc, const char *argv[]) {
  if (argc < 5) {
    Die("USAGE: bam_driver <input>.bam <child SM> <mother SM> <father SM>");
  }

  const string input = argv[1];
  const string child_sm = argv[2];
  const string mother_sm = argv[3];
  const string father_sm = argv[4];

  BamReader reader;
  reader.Open(input);
  if (!reader.IsOpen()) {
    Die("Input file could not be opened.");
  }
  TrioVector sites = BamUtility::ParseSites(reader, child_sm, mother_sm, father_sm);
  reader.Close();

  cout.precision(16);
  TrioModel params;
  ParameterEstimates *stats = EstimateParameters(params, sites);
  stats->IsLogLikelihoodIncreasing();
  stats->PrintMaxSequencingErrorRateEstimate();

  return 0;
}
