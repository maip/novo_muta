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
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/mip/novo_muta_infinite_sites_model/bamtools/lib -I/home/mip/novo_muta_infinite_sites_model/bamtools/include -lbamtools -I/home/mip/novo_muta_infinite_sites_model/bamtools/src -lbamtools-utils -I/usr/local/include -o bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc em_algorithm.cc sufficient_statistics.cc bam_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./bam_driver <output>.txt <input>.bam <child SM> <mother SM> <father SM>
 */
#include "sufficient_statistics.h"
#include "variant_visitor.h"


int main(int argc, const char *argv[]) {
/*  if (argc < 5) {
    Die("USAGE: bam_driver <input>.bam <child SM> <mother SM> <father SM>");
  }

  const string input = argv[1];
  const string child_sm = argv[2];
  const string mother_sm = argv[3];
  const string father_sm = argv[4];
  const int qual_cut = 13;  // May decide to pass cut values via command line.
  const int mapping_cut = 13;
  const double probability_cut = 0.0;  // Examines all probabilities, otherwise 0.1

  BamReader reader;
  reader.Open(input);
  if (!reader.IsOpen()) {
    Die("Input file could not be opened.");
  }

  BamAlignment al;
  TrioModel params;
  PileupEngine pileup;
  SamHeader header = reader.GetHeader();
  RefVector references = reader.GetReferenceData();
  VariantVisitor *v = new VariantVisitor(references, header, params, al,
                                         child_sm, mother_sm, father_sm,
                                         qual_cut, mapping_cut, probability_cut);
  
  pileup.AddVisitor(v);
  while (reader.GetNextAlignment(al)) {
    pileup.AddAlignment(al);
  }
  
  pileup.Flush();
  reader.Close();
*/
  // EM algorithm begins with initial E-Step.
  // Runs EM on a section created by a trio repeated 10x for every trio.
  //ReadDataVector data_vec = {{0, 0, 0, 4}, {0, 0, 4, 0}, {0, 1, 0, 3}};
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  int nan_count = 0;
  int one_count = 0;
  for (auto data_vec : trio_vec) {
    TrioModel params;
    TrioVector sites;
    for (int i = 0; i < 10; ++i) {
      sites.push_back(data_vec);
    }

    int sites_count = sites.size();  // 10
    SufficientStatistics stats(sites_count);
    stats.Update(params, sites);

    // M-Step.
    int count = 0;
    double maximized = stats.MaxSequencingErrorRate();
    double log_likelihood = stats.log_likelihood();
    cout.precision(16);

    //cout << "Performing the EM algorithm..." << endl;
    while (!Equal(params.sequencing_error_rate(), maximized) &&
        count < 50 && !stats.IsNan()) {  // Exits if converges or takes longer than 50 iteratons.
      // Sets new estimate.
      //cout << count << " ITERATION" << endl;
      params.set_sequencing_error_rate(maximized);
      stats.Clear();  // Sets all statistics except number of sites to 0.
      stats.Update(params, sites);  // Loops to E-Step.
      //stats.Print();

      // Sum of likelihood should increase and converge.
      if (stats.log_likelihood() < log_likelihood) {
        cout << "ERROR: Log likelihood is decreasing between iterations." << endl;
        PrintReadDataVector(params.read_dependent_data().read_data_vec);
        stats.Print();
      }

      if (stats.IsNan()) {
        nan_count++;
        cout << "ERROR: This site is nan." << endl;
        stats.Print();
      }

      maximized = stats.MaxSequencingErrorRate(); // Loops to M-Step.
      if (maximized == 1.0) {
        one_count++;
      }
      //cout << "~E:\t" << maximized << endl;
      log_likelihood = stats.log_likelihood();
      count++;
    }

    if (stats.IsNan()) {
      PrintReadDataVector(params.read_dependent_data().read_data_vec);
    }
    /*cout << "After " << count << " iterations of "
         << sites_count << " sites:" << endl
         << "^E:\t" << params.sequencing_error_rate() << endl << endl;*/
  
    stats.Clear();
  }
  cout << nan_count << endl;
  cout << one_count << endl;

  /*TrioVector sites = v->sites();
  int sites_count = sites.size();
  if (sites_count > 0) {
    SufficientStatistics stats(sites_count);
    stats.Update(params, sites);

    // M-Step.
    int count = 0;
    double maximized = stats.MaxSequencingErrorRate();
    double log_likelihood = stats.log_likelihood();
    cout.precision(16);
    while (!Equal(params.sequencing_error_rate(), maximized) &&
        count < 50) {  // Exits if converges or takes longer than 50 iteratons.
      params.set_sequencing_error_rate(maximized);  // Sets new estimate.
      stats.Clear();  // Sets all statistics except number of sites to 0.
      stats.Update(params, sites);  // Loops to E-Step.

      cout << "~E:\t" << params.sequencing_error_rate() << endl;
      stats.Print();

      // Sum of likelihood should increase and converge.
      if (stats.log_likelihood() < log_likelihood) {
        cout << "ERROR: Log likelihood is decreasing between iterations." << endl;
      }

      maximized = stats.MaxSequencingErrorRate(); // Loops to M-Step.
      log_likelihood = stats.log_likelihood();
      count++;
    }
    
    cout << "After " << count << " iterations of "
         << sites_count << " sites:" << endl
         << "^E:\t" << params.sequencing_error_rate() << endl;
  }*/

  return 0;
}
