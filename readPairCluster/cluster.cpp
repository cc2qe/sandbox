#include <vector>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <cstdlib>
#include <ios>

#include <unistd.h>

#include "simpleGraph.hpp"
#include "readPair.hpp"

using namespace std;

// input file format
// HWI-EAS87:7:30:1550:1386#0      chr1    3181533 97      chr10   102552618       145     76      37      76      37

typedef vector<string *> * tline;


// read all read pairs into a list
// 
list<readPair *> * readFile(string filename) {
  ifstream f(filename.c_str());
  list<readPair *> * newList = new list<readPair *>;
  string line;
  stringstream ss;
  string * datum;
  string seg;
  long linecount = 0;
  readPair *prev = NULL;

  while (!f.eof()) {
    linecount++;
    tline v = new vector<string *>;
    getline(f, line);
    ss << line;
    if (linecount % 10000 == 0)  cerr << linecount << endl; //<< line << endl;
    while (getline(ss,seg,'\t')) {
      datum = new string(seg);
      // filter out duplicates here
      v->push_back(datum);
    }
    
    if (v->size() > 0) { 
      readPair *rp = new readPair(*v);
      if (prev == NULL || !(prev->chrA == rp->chrA && prev->chrB == rp->chrB && prev->posA == rp->posA && prev->posB == rp->posB))
      newList->push_back(rp);
      prev = rp;
      
      // free the data
      
    }
    //	cout << rp->readID << rp->posA << endl;
    ss.clear();
  }
  return newList;
}

// load a list of read pairs from the list that are potential candidates to be clustered together
int getReadPairGroup (ifstream &f, long &startLoc, list<readPair* > ** dL, long DIST_CUTOFF) {
  // long prevLinePos, thisLinePos;
  // string *datum;
  // string line, seg;
  // bool inBounds = true;
  // stringstream ss;
  // long lastLoc, linecount;
  //  list<readPair *> * newList = new list<readPair *>;

  list<readPair *> * newList = new list<readPair *>;

  string line;
  stringstream ss;
  string * datum;
  string seg;
  long linecount = 0;
  long lastLoc = 0;
  readPair *prev = NULL;


  f.seekg(startLoc);

  while (!f.eof()) {
    ss.clear();
    lastLoc = f.tellg();
    linecount++;
     tline v = new vector<string *>;

     getline(f, line);
     ss << line;
     
     // cout << "current position: " << f.tellg() << "\n";

     // tokenize the line and push it into an array
     while (getline(ss, seg, '\t')) {
       datum = new string(seg);
       v->push_back(datum);
     }

     //     cout <<  v->size() << "\n";
     
    if (v->size() > 0) {
      // cout << "got a line!  \n";
      readPair *rp = new readPair(*v);
      
      // free memory
      for (vector<string *>::iterator it = v->begin(); it != v->end(); it++) {
	delete (*it);
      }
      delete v;

      if (prev == NULL || !(prev->chrA == rp->chrA && prev->chrB == rp->chrB && prev->posA == rp->posA && prev->posB == rp->posB)) {
  	// if it is within bounds, then add it to the list
  	if (prev == NULL || (prev->chrA == rp->chrA && prev->chrB == rp->chrB && ((rp->posA - prev->posA) < DIST_CUTOFF) )) {
  	  newList->push_back(rp);
  	  prev = rp;
  	}
  	else {
  	// else return the current list early, set the file pointer to just before the last line read
  	  startLoc = lastLoc;
  	  delete rp;
	  *dL = newList;
  	  return 1;
  	}
      }
    }
    

  }
  return 0; // newList;
}



void printUsage() {
  cerr << "readPairCluster" << endl;
  cerr << "-----------------------------" << endl;
  cerr << "usage:  readPairCluster -r readpairfile.txt [-d DIST_CUTOFF] [-q MAPQ_THRESH] [-s MIN_CLUST_SIZE] [-u]" << endl;
  cerr << "options: " << endl;
  cerr << "  -d DISTANCE_CUTOFF :: max distance between ends of reads that are clustered together (default 3kb)" << endl;
  cerr << "  -q MAPQ_THRESHOLD :: quality score cutoff for including reads (default 50)" << endl;
  cerr << "  -s MIN_CLUST_SIZE :: minimum size of clusters to output (default 3 read pairs) " << endl;
  cerr << "  -u :: only consider unique reads " << endl;
  cerr << endl;
}


// clusters a list of readPairs that are potentially clusterable
// and write out the potential clusters
//
//  int DIST_CUTOFF = 3000;
//  int MAPQ_THRESHOLD = 50;
//  int MIN_RESULT_CLUSTER_SIZE = 3;
//  bool ELIM_DUPLICATES = false;

//
int clusterList(list<readPair *> * datalist, long * clusterNum, int DIST_CUTOFF, int MAPQ_THRESHOLD, int MIN_RESULT_CLUSTER_SIZE, bool ELIM_DUPLICATES) {

  list<readPair *>::iterator first, last, next, r, s, rprev, sprev;

  long length = distance(datalist->begin(), datalist->end());
  long curCount = 0;

    first = datalist->begin();
    last = first;
    
    // start iterating through all pairs 
    while (first != datalist->end()) {
  
      for (next = first;
	   next != datalist->end() 
	   	&& abs( (*last)->posA - (*next)->posA ) <  DIST_CUTOFF  // why 
	   	 && (*last)->chrA == (*next)->chrA
	   	 && (*last)->chrB == (*next)->chrB;
	   next++) {
	last = next;
	curCount++;
      }
    
      if (distance(first, last) > 0) {
	simpleGraph<list<readPair *>::const_iterator> g = simpleGraph<list<readPair *>::const_iterator>();
	map< pair<long, long>, int > coordBox;
	rprev = first;

	for (r = first; r != last; r++) {

	  rprev = r;
	  coordBox[pair<long, long>( (*r)->posA, (*r)->posB)] = 1;

	  if (rprev == r || (*r)->posA != (*rprev)->posA || (*r)->posB != (*rprev)->posB) {
	  sprev = r;
	  for ( s = r; s != next; s++) {

	    if (s != r && ((*r)->chrB.compare((*s)->chrB) == 0) 
		&& abs( (*r)->posB - (*s)->posB ) < DIST_CUTOFF
		&& abs( (*r)->posA - (*s)->posA ) < DIST_CUTOFF
		&& (*r)->chrB == (*s)->chrB
		&& ( (*s)->posA != (*sprev)->posA || (*s)->posB != (*sprev)->posB ) 
		&& ( (*r)->posA != (*s)->posA || (*r)->posB != (*s)->posB )
		&& (*r)->scoreA > MAPQ_THRESHOLD 
		&& (*r)->scoreB > MAPQ_THRESHOLD 
		&& (*s)->scoreA > MAPQ_THRESHOLD
		&& (*s)->scoreB > MAPQ_THRESHOLD ) {

	      if (coordBox.find(pair<long, long>( (*s)->posA, (*s)->posB)) == coordBox.end()) {
			g.addEdge(r, s);
	      }

	      if (ELIM_DUPLICATES) {
			coordBox[pair<long, long>( (*s)->posA, (*s)->posB)] = 1;
	      }
	    }
	    
	    sprev = s;
	  }
	  }
	}
	
      g.connectedComponents(MIN_RESULT_CLUSTER_SIZE, clusterNum);
//      cout << distance(first, last) << "\t" << (*first)->chrA << "\t" << (*last)->chrA << "\t" << (*first)->chrB << "\t" << (*last)->chrB << " ---------------- " << endl << endl;

      }
        
    first = next;
    last = first;
    
    }
  }



int main(int argc, char **argv) {

  // get command-line options
  char *cvalue = NULL;
  int c;
  string readFileName;

  int DIST_CUTOFF = 3000;
  int MAPQ_THRESHOLD = 50;
  int MIN_RESULT_CLUSTER_SIZE = 3;
  bool ELIM_DUPLICATES = false;
 
  while ((c = getopt(argc, argv,"r:d:q:s:u")) != -1) {
    switch (c)
      {
      case 'r':
	readFileName = optarg; 
	//cout << readFileName << endl;
	break;
      case 'd':
	DIST_CUTOFF = atoi(optarg);
	break;
      case 'q':
	MAPQ_THRESHOLD = atoi(optarg);
	break;
      case 's':
	MIN_RESULT_CLUSTER_SIZE = atoi(optarg);
	break;
      case 'u':
	ELIM_DUPLICATES = true;
	break;
      case '?':
	if (optopt == 'r') {
	  printUsage();
	  cerr << "Option -r requires an argument (the filename containing the reads)" << endl;
	  }
      default:
	printUsage();
	cerr << "readPairCluster" << endl;
	abort();

      }
  }

  if (readFileName.size() == 0) {
      printUsage();
      cerr << "No read pair file specified!" << endl;
      abort();
  }
  
  cerr << "Only including reads better than MAPQ threshold of " << MAPQ_THRESHOLD << endl;
  cerr << "Using distance cutoff for clustering of " << DIST_CUTOFF << endl;
  cerr << "Min result cluster size " << MIN_RESULT_CLUSTER_SIZE << endl;


  long startLine = 0;
  ifstream f(readFileName.c_str());

  list<readPair *> *dL;
  long clusterNum = 0;


  while (getReadPairGroup(f, startLine, &dL, DIST_CUTOFF)) {
    clusterList(dL, &clusterNum, DIST_CUTOFF, MAPQ_THRESHOLD, MIN_RESULT_CLUSTER_SIZE, ELIM_DUPLICATES);
    //    cerr << startLine << "\t" << dL->size() << "\n";

    // need to iterate through this list and delete all contents??
    for (list<readPair *>::iterator it = dL->begin(); it != dL->end(); it++) {
      delete *it;
    }
    delete dL;

  }


  //  list<readPair *> * datalist = readFile(readFileName);    
  //  cerr << "Finished loading data file...." << endl;
 
  //  clusterList(datalist, DIST_CUTOFF, MAPQ_THRESHOLD, MIN_RESULT_CLUSTER_SIZE, ELIM_DUPLICATES);
 

}
