/**
 * @file indiv_driver.cc
 * @author Melissa Ip
 *
 * This allows the user to pass in an individual site that is converted to a
 * ReadDataVector. Using the trio model, the probability of mutation is
 * calculated and printed.
 */
#include <sstream>

#include "multinomial_trio_model.h"
#include "trio_model.h"


int main() {
  string child_str;
  string mother_str;
  string father_str;
  uint16_t A = 0;
  uint16_t C = 0;
  uint16_t G = 0;
  uint16_t T = 0;

  // char is_dirichlet;
  // cout << "Select a for multinomial trio model or b for Dirichlet-multinomial trio model." << endl;
  // cin.get(is_dirichlet);
  // cin.ignore(20, '\n');  // Flush buffer.

  // if (is_dirichlet == 'b') {
  //   TrioModel params;
  // } else {
  //   MultinomialTrioModel params;
  // }

  TrioModel params;
  
  char is_default;
  cout << "Use default trio model parameters? y/n" << endl;
  cin.get(is_default);
  cin.ignore(20, '\n');  // Flush buffer.

  if (is_default == 'n') {  // All other values uses default parameters.
    string rate_str;
    double rate = 0.0;
    cout << endl << "Enter population mutation rate (theta):" << endl;
    getline(cin, rate_str);
    stringstream rate_ss1(rate_str);
    rate_ss1 >> rate;
    params.set_population_mutation_rate(rate);

    cout << endl << "Enter germline mutation rate:" << endl;
    getline(cin, rate_str);
    stringstream rate_ss2(rate_str);
    rate_ss2 >> rate;
    params.set_germline_mutation_rate(rate);

    cout << endl << "Enter somatic mutation rate:" << endl;
    getline(cin, rate_str);
    stringstream rate_ss3(rate_str);
    rate_ss3 >> rate;
    params.set_somatic_mutation_rate(rate);

    cout << endl << "Enter sequencing error rate:" << endl;
    getline(cin, rate_str);
    stringstream rate_ss4(rate_str);
    rate_ss4 >> rate;
    params.set_sequencing_error_rate(rate);
  }

  while (true) {
    cout << endl << "Enter the child sequencing read <A C G T>. Example: 0 1 2 3" << endl;
    getline(cin, child_str);
    stringstream child(child_str);
    child >> A;
    child >> C;
    child >> G;
    child >> T;
    ReadData child_read = {A, C, G, T};

    cout << endl << "Enter the mother sequencing read <A,C,G,T>. Example: 0,1,2,3" << endl;
    getline(cin, mother_str);
    stringstream mother(mother_str);
    mother >> A;
    mother >> C;
    mother >> G;
    mother >> T;
    ReadData mother_read = {A, C, G, T};

    cout << endl << "Enter the father sequencing read <A,C,G,T>. Example: 0,1,2,3" << endl;
    getline(cin, father_str);
    stringstream father(father_str);
    father >> A;
    father >> C;
    father >> G;
    father >> T;
    ReadData father_read = {A, C, G, T};

    ReadDataVector data = {child_read, mother_read, father_read};
    cout << endl << "P(Mut):\t" << params.MutationProbability(data) << endl << endl;
  }

  return 0;
}
