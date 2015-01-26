/**
 * @file parse_bam_header.cc
 * @author Melissa Ip
 *
 * This reads the alignments of a bam file and parses read group data marked
 * with a @RG tag and outputs it as a tap separated text file with each read
 * group on a new line. The input for this file will be the outputed header
 * data from the bam file. This step is used to create a read group file in
 * order to index a bam region slice and is part of the input preparation for
 * bam_driver.cc. Repeat for each family member bam file.
 *
 * To merge all three read group files together:
 * cat <child_rg>.txt <mother_rg>.txt <father_rg>.txt > <rg>.txt
 *
 * To compile on Herschel:
 * c++ -std=c++11 -L/usr/local/lib -I/usr/local/include -o parse_bam_header utility.cc parse_bam_header.cc
 *
 * To run this file, provide the following command line inputs:
 * ./parse_bam_header <input>.txt <output>.txt
 */
#include <fstream>
#include <sstream>

#include "utility.h"


int main(int argc, const char *argv[]) {
  if (argc < 3) {
    Die("USAGE: parse_bam_header <input>.txt <output>.txt");
  }

  const string file_name = argv[1];
  ifstream fin(file_name);
  if (!fin.is_open() || 0 != fin.fail()) {
    Die("Input file cannot be read.");
  }

  const string fout_name = argv[2];
  ofstream fout(fout_name);

  string line;
  string tag;
  while (getline(fin, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    stringstream str(line);
    str >> tag;
    if (tag.compare("@RG") == 0) {
      fout << line << endl;
    }
  }
  fin.close();
  fout.close();

  return 0;
}
