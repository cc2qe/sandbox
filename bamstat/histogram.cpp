#include "histogram.hpp"

histogram::histogram(float minValue, float maxValue, int numBins) {
    this->binWidth = (maxValue - minValue) / float(numBins); 
    this->underBin = 0;
    this->overBin = 0;
    this->numBins = numBins;
    this->bins.resize(numBins);
    
    this->minValue = minValue;
    this->maxValue = maxValue;
    
    for (int i = 0; i < this->numBins; i++) {
        this->bins[i] = 0;
    }
    
}

int histogram::count(float dataval) {
    if (dataval < minValue) { 
        underBin++;
    } else if (dataval > maxValue) { 
        overBin++;
    } else {
        int binNum = int((dataval - this->minValue) / this->binWidth);
        bins[binNum]++;
    }
    
    return 0;
}

int histogram::writeFile(std::string filename) {
    std::ofstream histfile(filename.c_str());

    std::cout << filename;
    
    histfile << "# histogram" << std::endl;
    histfile << "# min of bin range " << this->minValue << std::endl;
    histfile << "# max of bin range " << this->maxValue << std::endl;
    
    for (int i = 0; i < this->numBins; i++) {
      histfile << i << "\t" << bins[i] << "\t" << (minValue + i*this->binWidth) << std::endl; // "->" << (minValue + (i+1)*this->binWidth) << "\t" << std::endl;
    }
    
    histfile.close();
}

histogram::~histogram() {

}
