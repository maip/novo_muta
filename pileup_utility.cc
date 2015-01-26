/**
 * @file pileup_utility.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the functions needed to parse and
 * manipulate pileup data.
 *
 * See top of pileup_utility.h for a complete description.
 */
#include "pileup_utility.h"
 

/**
 * Removes initial sequences from the ifstream that are not necessary, which
 * contain a N reference (a site that is not matched to any base). Assumes that
 * all three pileup files contain the same number of invalid sequences. Returns
 * the first line where N is not the reference.
 *
 * @param  f Input stream.
 * @return   First line where N is not the reference.
 */
string PileupUtility::TrimHeader(ifstream &f) {
  string line;
  while (getline(f, line)) {
    line = GetSequence(line);  // Trims newline.
    if (!line.empty()) {
      return line;
    }
  }
  return "";  // ERROR: There were only invalid N sequences.
}

/**
 * Trims the newline fron the end of the line, and returns the line if it is a
 * valid sequence that does not contain a N reference.
 *
 * @param  line Line from pileup file.
 * @return      Line without newline and has valid nucleotide reference.
 */
string PileupUtility::GetSequence(string &line) {
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  
  line.erase(remove(line.begin(), line.end(), '\n'), line.end());
  stringstream str(line);
  str >> sequence;
  str >> position;
  str >> ref_nucleotide;

  if (ref_nucleotide != 'N') {
    return line;
  } else {
    return "";
  }
}

/**
 * Parses pileup data into ReadData. Periods and commas match the
 * reference nucleotide.
 *
 * @param  line Read from a pileup file representing a single site sequence.
 * @return      ReadData.
 */
ReadData PileupUtility::GetReadData(const string &line) {
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  int num_aligned_reads = 0;
  string bases;

  stringstream str(line);
  str >> sequence;
  str >> position;
  // str.seekg(10, ios::beg);
  str >> ref_nucleotide;
  str >> num_aligned_reads;
  str >> bases;

  uint16_t matches = (count(bases.begin(), bases.end(), '.') +
                      count(bases.begin(), bases.end(), ','));
  uint16_t A = (count(bases.begin(), bases.end(), 'A') + 
                count(bases.begin(), bases.end(), 'a'));
  uint16_t C = (count(bases.begin(), bases.end(), 'C') + 
                count(bases.begin(), bases.end(), 'c'));
  uint16_t G = (count(bases.begin(), bases.end(), 'G') + 
                count(bases.begin(), bases.end(), 'g'));
  uint16_t T = (count(bases.begin(), bases.end(), 'T') + 
                count(bases.begin(), bases.end(), 't'));

  ReadData data;
  if (ref_nucleotide == 'A') {
    data = {matches, C, G, T};  
  } else if (ref_nucleotide == 'C') {
    data = {A, matches, G, T}; 
  } else if (ref_nucleotide == 'G') {
    data = {A, C, matches, T}; 
  } else if (ref_nucleotide == 'T') {
    data = {A, C, G, matches};
  }
  return data;
}

/**
 * Returns the ReadDataVector from the parsed pileup data.
 *
 * @param  child_line  Line from the child pileup.
 * @param  mother_line Line from the mother pileup.
 * @param  father_line Line from the father pileup.
 * @return             ReadDataVector from the parsed pileup data.
 */
ReadDataVector PileupUtility::GetReadDataVector(const string &child_line,
                                                const string &mother_line,
                                                const string &father_line) {
  ReadDataVector data_vec = {GetReadData(child_line),
                             GetReadData(mother_line),
                             GetReadData(father_line)};
  return data_vec;
}

/**
 * Opens and parses all pileup files. All valid sequences are converted to
 * ReadDataVector.
 *
 * @param  child  Chile pileup input stream.
 * @param  mother Mother pileup input stream.
 * @param  father Father pileup input stream.
 * @return        List of probabilities from pileup data.
 */
TrioVector PileupUtility::ParseSites(ifstream &child,
                                     ifstream &mother,
                                     ifstream &father) {
  // Removes N sequences.
  string child_line = TrimHeader(child);
  string mother_line = TrimHeader(mother);
  string father_line = TrimHeader(father);
  if (child_line.empty() || mother_line.empty() || father_line.empty()) {
    Die("Pileup file does not contain valid sequences (no N reference).");
  }

  TrioVector sites;
  ReadDataVector data_vec = GetReadDataVector(child_line, mother_line, father_line);
  sites.push_back(data_vec);  // First valid line.

  while (getline(child, child_line)) {
    getline(mother, mother_line);
    getline(father, father_line);
    data_vec = GetReadDataVector(child_line, mother_line, father_line);
    sites.push_back(data_vec);
  }

  return sites;
}

/**
 * Returns a list of probabilities of each site in pileup data at its sequence
 * position if it is greater than the threshold value.
 *
 * @param  sites List of sites from pileup data.
 * @return       List of probabilities from pileup data.
 */
vector<double> PileupUtility::GetProbability(const TrioVector &sites) {
  TrioModel params;
  vector<double> probabilities;
  for (auto site : sites) {
    double probability = params.MutationProbability(site);
    // if (probability >= kThreshold) {
      probabilities.push_back(probability);
    // }
  }
  return probabilities;
}

/**
 * The output file is tab separated and each sequence is on a new line.
 * The first column represents the sequence position and the second column
 * represents the probability at that sequence.
 *
 * @param  file_name     Output file name.
 * @param  probabilities List of probabilities from pileup data.
 */
void PileupUtility::WriteProbability(const string &file_name,
                                     const vector<double> &probabilities) {
  ofstream fout(file_name);
  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();
}
