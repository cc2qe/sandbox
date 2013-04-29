#include "discbins.hpp"


DiscBinFreqTable::DiscBinFreqTable(long min, long max) {
	// min to max are stored in counts[1] to counts[max - min + 1]
	// values less than min go in counts[0]
	// values more than min go into counts[max-min+2]
	
	this->min = min;
	this->max = max;
	this->numCounts = (max - min) + 3;
	this->totalValues = 0;
	this->counts = new long[this->numCounts];
	
	// initialize counts
	for (long i = 0; i < this->numCounts; i++) {
		this->counts[i] = 0;
	} 
}

long DiscBinFreqTable::count(long value) {
	this->totalValues++;
	
	if (value < this->min) {
		this->counts[0] += 1;
	} 
	else if (value > this->max) { 
		this->counts[this->numCounts + 1] += 1;
	} 
	else {
		this->counts[(value - min) + 1] += 1;
	}
	
}

long DiscBinFreqTable::getMedian() {
	long countSoFar = 0;
	long halfCount = this->totalValues / 2;
	long lastBin = 0;
	long i = 0;
	for (i = 0; countSoFar < halfCount; i++) {
		countSoFar += this->counts[i];
		if (this->counts[i] > 0) lastBin = i;
	}
	
	if (countSoFar == halfCount) { 
	  return (i-2) + this->min;
	  } else {
   	  return (lastBin-2) + this->min;	
	}
}

DiscBinFreqTable::~DiscBinFreqTable() {
	delete this->counts;
}
