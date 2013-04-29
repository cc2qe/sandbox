#include "readPair.hpp" 

readPair::readPair(vector<string*> &linedata) {
  // check that the line has the right number of values
  if (linedata.size() != 15) {
    throw "Wrong number of values in read pair line"; 
  }
  
  readID = *(linedata[0]);
  chrA = *(linedata[1]);
  chrB = *(linedata[4]);
  posA = atol((*linedata[2]).c_str());
  posB = atol((*linedata[5]).c_str());
  scoreA = atoi((*linedata[8]).c_str());
  scoreB = atoi((*linedata[10]).c_str());  
  lenA = atoi((*linedata[7]).c_str());
  lenB = atoi((*linedata[9]).c_str());  
  flagA = atoi((*linedata[3]).c_str());
  flagB = atoi((*linedata[6]).c_str());  
  
  seqA = *(linedata[11]);
  qualA = *(linedata[12]);

  seqB = *(linedata[13]);
  qualB = *(linedata[14]);

}

// readPair::readPair(bam1_t *a, bam1_t *b) {
//  readID = "readID";
//  chrA = 
//}

void readPair::print() {
  cout << this->readID << "\t";
  cout << this->chrA << "\t";
  cout << this->chrB << "\t";
  cout << this->posA << "\t";
  cout << this->posB << "\t";
  cout << this->scoreA << "\t";
  cout << this->scoreB << "\t";
  cout << this->lenA << "\t";
  cout << this->lenB << "\t";
  cout << this->flagA << "\t";
  cout << this->flagB << "\t";

  cout << this->seqA << "\t";
  cout << this->seqB << "\t";

  cout << this->qualA << "\t";
  cout << this->qualB << "\t";

  cout << endl;
  
}

readPair::~readPair() {
 
}
