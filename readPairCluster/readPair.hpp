#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "sam.h"

using namespace std;

class readPair {
public:
  readPair(vector<string*> &linedata);

  // make a readpair object from a pair of BAM records
  readPair(bam1_t *a, bam1_t *b);

  ~readPair();

  void print();

  // to keep things simple, all of these are left public
  string readID;
  string chrA;
  string chrB;
  long posA;
  long posB;
  int scoreA;
  int scoreB;
  int lenA;
  int lenB;
  int flagA;
  int flagB;

  string seqA;
  string seqB;
  
  string qualA;
  string qualB;


};
