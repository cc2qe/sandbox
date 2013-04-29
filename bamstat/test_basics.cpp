#include <iostream>
#include "discbins.hpp"

int main() {

	DiscBinFreqTable * db = new DiscBinFreqTable(9, 20);
	
	
	long values[16] = {2, 2, 8, 8, 4, 8, 12, 12, 8, 9, 10, 10, 10, 10, 11, 14};

	for (int i = 0; i < 16; i++) {
		db->count(values[i]);
	}
	
	std::cout << "The median is " << db->getMedian();
	
	delete db;
 
 	return 0;
}
